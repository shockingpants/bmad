def dec(fn,offset=0,sep=', '):
	##{{{
	"""
	Stores and retrive data of interest
	index is stored starting from 0. 0 usually points to the cap ACE.
	"""
	a=fn()
	assert type(a)==list
	def new_fn(offset=0,sep=' ,'):
		dat=sep.join(map(lambda x: str(x+offset),a))
		return dat
	return new_fn
	##}}}

@dec
def sec_site():
	##{{{
	"""
	Within 4A
	"""
	dat=[1,2,4,25,26,76,77,79,80,83,85]
	return dat
	##}}}

###############################
# Usage: 
sec_site()
sec_site(sep=', ')
sec_site(offset=4)
###############################

