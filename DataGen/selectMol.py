import numpy 

class EFGSample:
    """
    This is a EFGSample class. 
    The EFGSample class includes the basic sample functions.
    This class is used to do the sampling from large dataset based on EFGs.
    
    Parameters
    ----------
    input_dict : dictionary
        input_dict is the input EFG directory,
        with key as EFGs name and value as list of indexes of molecules containing this EFG.
    cutoff: int, optional
        the number of molecules for each EFG will be sampled in each iteration, by default 1
    numbers : int, optional
        the number of molecules will be sampled, by default 99000. Either numbers or iterations should be given.
    iterations : int, optional
        the number of iterations will be looped for whole EFGs, by default None
    sort : bool, optional
        sort the input_dict based on length of its values, by default True
    previous_list : list, optional
        list of molecules will be not processed, by default None
    """
    def __init__(self, input_dict, cutoff=1, numbers=None, iterations=None, sort=True, previous_list=None):
        super(EFGSample, self).__init__()
        self.input_dict = input_dict
        self.cutoff = cutoff
        self.numbers = numbers
        self.iterations = iterations
        self.sort = sort
        self.previous_list = previous_list

        # update input_dict by removing all molecules in previous_list #
        if not self.previous_list:
            self.previous_list = []
        else:
            for index in self.previous_list:
                for value in self.input_dict.values():
                    if index in value:
                        value.remove(index)

        # sort input_dict based on length of values 
        if self.sort:
            self.sorted_input_dict = sorted(self.input_dict.items(), key=lambda item: len(item[1]), reverse=True)
        else:
            self.sorted_input_dict = [[key,value] for key, value in self.input_dict.items()]

    @property
    def sampled_index(self):
        """
        sampled_index is a list used to save sampled molecule indexes

        """
        self._sampled_index = []
        # self._sampled_index_before is used to save molecules sampled in all previous iterations #
        self._sampled_index_before = []
        self._sampled_EFGs_dic = {key[0]:[] for key in self.sorted_input_dict}
        if self.iterations:
            self._sample_iters()
            return self._sampled_index
        if self.numbers:
            self._sample_numbers()
            return self._sampled_index
            
    @property
    def sampled_EFGs_dic(self):
        """
        sampled_EFGs_dic is a dictionary used to save sampled molecule indexes based on EFGs
        sampled_index should be first conducted before calling sampled_EFGs_dic
        """
        try:
            return self._sampled_EFGs_dic
        except:
            raise NotImplementedError('sampled_index should be first conducted')

    def _sample_basic(self, index_list, name):
        """
        basic sampling function used in iteration sampling and numbers sampling
        
        Parameters
        ----------
        index_list : list
            list of indexes belonging to the given EFG name 
        name : str
            EFG name
        
        Returns
        -------
        append_list: list
            list of sampled indexes for the given index_list
        """
        append_list = []
        index_list = [i for i in index_list if (i not in self._sampled_index_before) and (i not in self.previous_list)]
        if len(index_list) >= self.cutoff:
            # if len(index_list) is larger than self.cutoff, we only sampled first self.cutoff molecules #
            append_list = [ i for i in index_list[0:self.cutoff]]
            self._sampled_index.extend(append_list)
            self._sampled_EFGs_dic[name].extend(append_list)
        elif len(index_list) != 0:
            # if len(index_list) is less than self.cutoff, we sampled all molecules #
            append_list = [ i for i in index_list]
            self._sampled_index.extend(append_list)
            self._sampled_EFGs_dic[name].extend(append_list)
        return append_list

    
    def _sample_iters(self):
        """
        iterations sampling 
        for each iteration, we sampled molecules no larger than cutoff for each EFG
        """
        for i in range(self.iterations):
            print("iters:" + str(i + 1))
            print("number of sampled index:" + str(len(self._sampled_index)))
            for EF in self.sorted_input_dict:
                index_list = EF[1]
                name = EF[0]
                _ = self._sample_basic(index_list, name)
            # update self._sampled_index by removing redudent molecules (keep order) after each iteration #
            self._sampled_index = list(dict.fromkeys(self._sampled_index))
            self._sampled_index_before = [t for t in self._sampled_index]

    
    def _sample_numbers(self):
        """
        numbers sampling 
        for each iteration, we sampled molecules no larger than cutoff for each EFG
        after each EFG sampling, we check the numbers of sampled molecules, 
        if it's equal or larger than self.number, stop sampling
        """
        i = 0 
        # continue until the number of sampled indexes is equal or larger the given self.numbers #
        while (len(self._sampled_index) < self.numbers):
            i += 1
            print("iters:" + str(i))
            print("number of sampled index:" + str(len(self._sampled_index)))

            for EF in self.sorted_input_dict:
                index_list = EF[1]
                name = EF[0]
                append_list = self._sample_basic(index_list, name)
                # update _sampled_index by removing redudent molecules (keep order) after each adding #
                self._sampled_index = list(dict.fromkeys(self._sampled_index))
                # check number #
                if len(self._sampled_index) >= self.numbers:
                    for index in append_list:
                        if index not in self._sampled_index[0:self.numbers]:
                            self._sampled_EFGs_dic[name].remove(index)  
                    self._sampled_index = self._sampled_index[0:self.numbers]
                    break    
            self._sampled_index_before = [t for t in self._sampled_index]

class FreqEFGSample(EFGSample):
    """
    This is a EFGSample class based on frequency.
    
    Parameters
    ----------
    input_dict : dictionary
        input_dict is the input EFG directory, with key as EFGs name and value as list of indexes of molecules containing this EFG.
    frequency_list: list
        frequency_list is the list of molecule frequencies
        the order of frequency numbers is determined by molecule index
        for example, the frequency of molecule with index 1 is frequency_list[1] 
    cutoff: int, optional
        the number of molecules for each EFG will be sampled in each iteration, by default 1
    numbers : int, optional
        the number of molecules will be sampled, by default 99000. Either numbers or iterations should be given.
    iterations : int, optional
        the number of iterations will be looped for whole EFGs, by default None
    sort : bool, optional
        sort the input_dict based on length of its values, by default True
    previous_list : list, optional
        list of molecules will be not processed, by default None

    """
    def __init__(self, input_dict, frequency_list, cutoff=1, numbers=None, iterations=None, sort=True, previous_list=None):
        super().__init__(input_dict, cutoff=cutoff, numbers=numbers, iterations=iterations, sort=sort, previous_list=previous_list)
        # sort index list for each key using index frequency #
        self.frequency_list = frequency_list
        sorted_input_dict_freq = []
        for item in self.sorted_input_dict:
            index_list = item[1]
            index_list_sorted = sorted(index_list, key=lambda index: self.frequency_list[index], reverse=True)
            sorted_input_dict_freq.append([item[0], index_list_sorted])
        self.sorted_input_dict = sorted_input_dict_freq


class SizeEFGSample(EFGSample):
    """
    This is a EFGSample class based on molecule size.
    
    Parameters
    ----------
    input_dict : dictionary
        input_dict is the input EFG directory, with key as EFGs name and value as list of indexes of molecules containing this EFG.
    size_list: list
        size_list is the list of molecule sizes
        the order of molecule sizes is determined by molecule index
        for example, the size of molecule with index 1 is size_list[1] 
    cutoff: int, optional
        the number of molecules for each EFG will be sampled in each iteration, by default 1
    numbers : int, optional
        the number of molecules will be sampled, by default 99000. Either numbers or iterations should be given.
    iterations : int, optional
        the number of iterations will be looped for whole EFGs, by default None
    sort : bool, optional
        sort the input_dict based on length of its values, by default True
    previous_list : list, optional
        list of molecules will be not processed, by default None

    """
    def __init__(self, input_dict, size_list, cutoff=1, numbers=None, iterations=None, sort=True, previous_list=None):
        super().__init__(input_dict, cutoff=cutoff, numbers=numbers, iterations=iterations, sort=sort, previous_list=previous_list)
        ### sort indexs for each key using size_list
        self.size_list = size_list
        self.size = max(self.size_list)

    def _sample_basic(self, index_list, name):
        """
        update basic sampling function by considering molecule size
        
        Parameters
        ----------
        index_list : list
            list of indexes belonging to the given EFG name 
        name : str
            EFG name
        
        Returns
        -------
        append_list: list
            list of sampled indexes for the given index_list
        """
        index_list = [i for i in index_list if (i not in self._sampled_index_before) and (i not in self.previous_list)]
        append_list = []
        # sampled molecules spread in molecule size #
        if len(index_list) >= self.size:
            size_dic = {idx: self.size_list[idx] for idx in index_list}
            size_dic_sorted = sorted(size_dic.items(), key=lambda item: item[1])
            min_size = min([i[1] for i in size_dic_sorted])
            max_size = max([i[1] for i in size_dic_sorted])
            while(len(append_list) < self.size):
                for s in range(min_size, max_size + 1):
                    for ss in size_dic_sorted:
                        if ss[1] == s and ss[0] not in append_list:                                  
                            append_list.append(ss[0])
                            break
            self._sampled_index.extend(append_list)
            self._sampled_EFGs_dic[name].extend(append_list)
        elif len(index_list) != 0:
            append_list = [ i for i in index_list]
            self._sampled_index.extend(append_list)
            self._sampled_EFGs_dic[name].extend(append_list)
        
        return append_list