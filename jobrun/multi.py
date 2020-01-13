from multiprocessing import Pool
import multiprocessing, logging
mpl = multiprocessing.log_to_stderr()
mpl.setLevel(logging.INFO)
import traceback
import pdb
class Multi:
    def multi_run_wrapper(twargs):
        fun=twargs['func']
        data=twargs['data']
        params=twargs['params']
        results=fun(**data,**params)
        return(results)   
