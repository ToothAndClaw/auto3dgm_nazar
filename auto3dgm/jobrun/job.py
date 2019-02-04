class Job(object):
	"""Data structure encapsulating data and parameters for JobRun task.

	Example dict version of Job data
	{
    ‘data’: 
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

	"""

	def __init__(self, job_dict={}, data={}, params={}, func=None):
		self.data = {}
		self.params = {}
		self.func = None

		if job_dict:
			self.import_job_dict(job_dict)
		elif data or params or func:
			self.import_args(data, params, func)

	def import_job_dict(self, job_dict):
		if (job_dict and isinstance(job_dict, dict)):
			if 'data' in job_dict and self.__validate_data(job_dict['data']):
				self.data = job_dict['data']
			if 'params' in job_dict and self.__validate_params(job_dict['params']):
				self.params = job_dict['params']
			if 'func' in job_dict and self.__validate_func(job_dict['func']):
				self.func = job_dict['func']

	def import_args(data={}, params={}, func={}):
		if data and __validate_data(data):
            self.data = data
	        
        if params and __validate_params(params):
            self.params = params
        
        if func and __validate_func(func):
            self.func = func
        
	def as_dict(self):
		"""Returns job data structure as dict"""
		if (__validate_data(self.data) 
			and __validate_params(self.params) 
			and __validate_func(self.func)):
		return {
			'data': self.data,
			'params': self.params,
			'func': self.func
		}

	def validate(self):
    	"""Check all components and return true if all validate"""
    	if (self.data and self.__validate_data(self.data) 
    		and self.params and self.__validate_params(self.params) 
    		and self.func and self.__vlaidate_func(self.func)):
    		return True

	def __validate_data(self, data):
		"""data must be dict, every element must be dict with >=1 element"""
        if (not data 
        	or not isinstance(data, dict) 
        	or not len(data)
        	or not self.__validate_data_items(data.values())):
            self.__validation_error(type='data', data)
        return True

    def __validate_data_items(self, items):
    	for x in items:
            if not isinstance(x, dict) or not len(x):
                self.__validation_error(type='data_item', x)
        return True

    def __validate_params(self, params):
        """Params must be dict with at least one value"""
        if not params or not isinstance(params, dict) or not len(params):
            self.__validation_error(type='params', params)
        return True

    def __validate_func(self, func):
        """Func must be callable"""
        if not func or not callable(func):
            self.__validation_error(type='func', func)
        return True  

    def __validation_error(self, error_type, var):
    	allowed_types = ['data', 'data_item', 'params', 'func']
    	if error_type not in allowed_types:
    		raise ValueError('Unexpected error type ' + str(error_type))
    	else:
    		raise ValueError('Unexpected value' + str(var) + 'for type ' + str(error_type))
