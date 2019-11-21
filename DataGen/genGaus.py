import os

def file_seperate(infile_list, datadir, outdir, initial_index=0):
    """ 
    Seperate conformations in *_confs.sdf 
    infile_list: a list of file name
    datadir: directory for input files
    outdir: directory for output seperate files
    initial_index: seperated files will be named as initial_index + number of conformations have been seperated . sdf
    return: seperated conformation files have been saved in outdir; return the (initial_index + the number of conformations have been seperated)
    """
    i = initial_index
    for name in infile_list:
        w = open(os.path.join(outdir,str(i + 1) + ".sdf"),"w")
        if name in os.listdir(datadir):
            infile = os.path.join(datadir, name)
        for line in open(infile):
            if line[0:4] == '$$$$':
                w.write(line)
                i += 1
                w.close()
                w = open(os.path.join(outdir,str(i+1) + '.sdf'), 'w')
            else:
                w.write(line)
        w.close()
    return i
    
def file_transfer(index_list, ftype, datadir):
    """ 
    Convert file into .com file using obabel 
    index_list: the list of all file name index that are needed to be converted 
    ftype: if cry --> local minimization, if others --> conformations
    datadir: the direcory for the sdf file, and all converted file will be in same directory
    notice: don't convert too many sdf files at the same time, it might cause problem of your computer
    """
    for i in index_list:
        if ftype == "cry":
            sdffn = str(i)+'_min.sdf'
        else:
            sdffn = str(i)+'.sdf'
        comfn =  str(i)+'.com'
        cmd = 'obabel -isdf ' + datadir + '/' + sdffn + ' -ocom -O ' + datadir + '/' + comfn
        os.system(cmd)

def gaussian_gen(index_list, header_file, datadir):
    """ 
    Add header for guassian optimization  
    header_file: the file with Gaussian calculation header information
    """
    header = header_file
    for i in index_list:
        comfn =  str(i)+'.com'
        outcomfn =  str(i)+'.opt.com'
        
        cmd1 = 'tail -n+4 ' + datadir + '/' + comfn + ' > ' + datadir + '/tmp' + str(i)
        cmd2 = 'cat ' + header + ' ' + datadir + '/tmp' + str(i) + ' > ' + datadir + '/' + outcomfn
        os.system(cmd1)
        os.system(cmd2)

