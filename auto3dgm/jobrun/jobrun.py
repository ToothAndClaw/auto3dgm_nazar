class JobRun(object):
    """Runs batch tasks in single core, multi core, or cluster environments 

   JobRun carries out tasks either locally or on a cluster by calling a 
   supplied function multiple times, with variable data objects and consistent 
   non-varying function parameters per call. Instances of this class must be 
   supplied with two sets of attributes: 1) local/cluster environment settings 
   (at the most basic level, this can just be a local/cluster switch flag), and 
   2) job parameters. Job parameters must include a function reference (which 
   itself should be able to receive named data variables and optionally further 
   named parameters), a dataset dictionary of data objects with keys indicating 
   individual job titles, and a dictionary of named parameters. Job parameters 
   should be able to be supplied either separately or as a unified data object: 

    {
    ‘dataset’: 
        {
        ‘analysis_1’: {‘mesh’: mesh1}, 
        ‘analysis_2’: {‘mesh’: mesh2}, 
        ‘analysis_3’: {‘mesh’: mesh3}
        }, 
    ‘params’: 
        {
        ‘point_number’: 200, 
        ‘subsample_method’: ‘GPR’
        }, 
    ‘func’: function_reference
    }

    Attributes:
        mode: Processing mode. String, can be 'single', 'multi', or 'cluster'
        dataset: Dict with analysis id keys and analysis data subdicts
        params: Dict for function args with arg name keys and arg values
        func: Function reference to be run for each analysis, using params
    """

    def __init__(self, analysis_dict={}, mode='', dataset={}, params={}, func=None, run=False):
        self.allowed_modes = ['single', 'multi']
        self.__mode = 'single'
        
        self.analysis_dict = {}
        self.dataset = {}
        self.params = {}
        self.func = None
        
        if mode and mode in self.allowed_modes:
            self.__mode = mode
        if analysis_dict:
            self.import_dict(analysis_dict)
        elif dataset and params and func:
            self.import_args(dataset, params, func)
        if run:
            self.execute_jobs()

    def import_dict(self, analysis_dict):
        if (analysis_dict and isinstance(analysis_dict, dict) and 
            self.__validate_dataset(dataset) and 
            self.__validate_params(params) and 
            self.__validate_func(func)
        ):
            self.analysis_dict = analysis_dict
            self.dataset = dataset
            self.params = params
            self.func = func
        else:
            msg = 'Problem validating analysis dict: {}'.format(analysis_dict)
            raise ValueError(msg)

    def import_args(self, mode='', dataset={}, params={}, func={}):
        if dataset and __validate_dataset(dataset):
            self.dataset = dataset
        else:
            raise ValueError('Problem validating dataset: {}'.format(dataset))

        if params and __validate_params(params):
            self.params = params
        else:
            raise ValueError('Problem validating params: {}'.format(params))

        if func and __validate_func(func):
            self.func = func
        else:
            raise ValueError('Problem validating func: {}'.format(func))

    def execute_jobs(self):
        if self.__mode and self.__mode in self.allowed_modes:
            if self.__mode == 'single':
                return run_single(self)
            elif self.__mode == 'multi':
                return run_multi(self)
            else:
                raise ValueError('Unexpected mode: {}'.format(self.__mode))
        else:
            raise ValueError('Current mode ({}) not an allowed mode: {}'.format(
                self.__mode, self.allowed_modes))

    def run_single(self):
        """Run jobs locally using a single core"""
        job_dict = {
            'output': {},
            'input': self.analysis_dict
        }

        results_dict = {}
        for k, v in self.analysis_dict.items():
            results_dict[k] = self.func(data=v, **self.params)

        job_dict['output']['results'] = results_dict
        return job_dict
        
    def run_multi(self):
        """Run jobs locally using mutliple cores"""

    def run_cluster(self):
        """Run jobs on a cluster environment (not yet implemented TODO TBA)"""

    # Private methods

    def __validate_dataset(self, dataset):
        """Dataset must be dict, every element must be dict with >=1 element"""
        success = True
        if dataset and isinstance(dataset, dict) and len(dataset):
            for x in dataset.values():
                if not isinstance(x, dict) or not len(x):
                    success = False
        else:
            success = False
        return success

    def __validate_params(self, params):
        """Params must be dict with at least one value"""
        if params and isinstance(params, dict) and len(params):
            return True
        else:
            return False

    def __validate_func(self, func):
        """Func must be callable"""
        if func and callable(func):
            return True
        else:
            return False

