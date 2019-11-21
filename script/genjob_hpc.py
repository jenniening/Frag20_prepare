import os

def gen_job(index_list, current_jobs, jobdir, datadir,  num_jobs=10, time=5, mem=10):
    """
    Generate jobs that can be submitted in HPC
    index_list: the list includes all indexes that we want to calculate
    current_jobs: the number of jobs that already be generated
    jobdir: the directory for saving job file
    datadir: the directory for the Gaussian input file and output file
    num_jobs: how many calculations that one job file includes, defaults to 10
    time: require hours from HPC for each job, defaults to 5
    mem: require memory from HPC for each job, defaults to 10
    notice: for large molecules, it shoule be better to contain only one caculation in each job, but for small molecules, contain more than one calculation in each job may be more efficient. Small time and mem can make it easier to get node to run the job. 
    """
    for n in range(int(len(index_list)/num_jobs) + 1):
        job = open(os.path.join(jobdir,"job_" + str(n) + ".pbs"), "w")
        job.write("#!/bin/bash\n#\n#SBATCH --job-name=G_" + initial_1 + "_" + str(n) + "\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=1\n#SBATCH --time=" + str(time) + "00:00\n#SBATCH --mem=" + str(mem) + "GB\n\nmodule purge\nmodule load gaussian/intel/g09e01\n\n")
        for i in range(num_jobs):
            i = n * num_jobs + i + current_jobs
            i = str(index_list[i])
            job.write("srun run-gaussian " + datadir + "/" + i + ".opt.com > " + datadir + "/" + i + ".opt.log 2>&1\n")
            job.close()


def sub_job(job_list, id1, id2):
    """
    Submit job 
    job_list: the list of jobs, this can be same as index_list when num_jobs=1, or can be a list of continous numbers 
    """
    for n in job_list[id1:id2]:
        print(n)
        os.system("sbatch job_" + str(n) + ".pbs")
    
def check_gaussian(index_list, datadir):
    """
    Check whether the gaussian calculation finishes normally 
    """
    error = open("Gaussian_error.csv", "w")
    for i in index_list:
        i = str(i)
        try: 
            infile = open(os.path.join(datadir, i + ".opt.log")
            line = [line for line in infile][-2]
            if line.split()[0:3] == ["Job","finishes","at:"]:
                continue
            else:
                error.write(i + "\n")
            infile.close()
        except:
            error.write(i + "\n")
    error.close()

