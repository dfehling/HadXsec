## Process is ttbar, qcd, or DATA
## Observable is taggedMass_all_pt, taggedMass_pt_400_500, taggedMass_pt_500_600, taggedMass_pt_600_700, taggedMass_pt_700_800,
##               taggedMass_pt_800


from theta_auto import *

def had_xs_model(filename, xsection):

	options = Options()
	options.set('global', 'debug', 'True')
	print options

	print filename

	ttbar_rate_unc = math.log(10.0)
	# qcd_rate_unc = math.log(1.5)

	## These two lines are for Barlow-Beeston lite
	ttbar_had_xs_model = build_model_from_rootfile(filename, include_mc_uncertainties=True)
	ttbar_had_xs_model.fill_histogram_zerobins(epsilon=None)

	# ttbar_had_xs_model = build_model_from_rootfile(filename, include_mc_uncertainties=False)
	# ttbar_had_xs_model.fill_histogram_zerobins()
	ttbar_had_xs_model.add_lognormal_uncertainty('ttbar_rate', ttbar_rate_unc, 'TTBAR')
	# ttbar_had_xs_model.add_lognormal_uncertainty('qcd_rate', qcd_rate_unc, 'QCD')

	# for p in ttbar_had_xs_model.distribution.get_parameters():
	# 	d = ttbar_had_xs_model.distribution.get_distribution(p)
		# if d['typ'] == 'gauss' and d['mean'] == 0.0 and d['width'] == 1.0:
		# 	ttbar_had_xs_model.distribution.set_distribution_parameters(p, range = [-5.0, 5.0])
		# if p=='ttbar_rate':
		# 	ttbar_had_xs_model.distribution.set_distribution(p, typ = "gauss", mean = 0.0, width = inf, range = [-inf,inf])
		# if p=='qcd_rate':
		# 	ttbar_had_xs_model.distribution.set_distribution(p, typ = "gauss", mean = 0.0, width = inf, range = [-inf,inf])

	# ttbar_had_xs_model.set_signal_processes('ttbar') ## ttbar signal
	signal_process_groups = {'': []} 
	# signal_process_groups = {'ttbar': ['ttbar']} ## ttbar signal

	parVals = mle(ttbar_had_xs_model, input = 'data', n=1, signal_process_groups = signal_process_groups)
	print parVals

	parameter_values = {}
	uncertainty_values = {}
	sf = {}
	for p in ttbar_had_xs_model.get_parameters([]):
		parameter_values[p] = parVals[''][p][0][0]
		uncertainty_values[p] = parVals[''][p][0][1]
	# 	# parameter_values[p] = parVals['ttbar'][p][0][0] ## ttbar signal
		# exponent = ttbar_rate_unc*parameter_values[p]#+math.pow(ttbar_rate_unc*uncertainty_values[p],2.0)/2.0
		# sf[p] = math.exp(exponent)
		# print p+" = "+str(sf[p])
	# parameter_values['beta_signal'] = parVals['ttbar']['beta_signal'][0][0] ## ttbar signal

	print parameter_values


	# res = ml_fit(model, signal_processes = signal_processes, **options)
	# for param in res[sp]:
	#            values[param] = res[sp][param][0][0]
	test = {}
	for proc in ttbar_had_xs_model.get_processes('taggedMass_2btag_000-040tau'):
		test[proc] = ttbar_had_xs_model.get_coeff('taggedMass_2btag_000-040tau', proc).get_value(parameter_values)
		if proc == "TTBAR":
			print test[proc]

	for proc in ttbar_had_xs_model.get_processes('taggedMass_1btag_000-040tau'):
		test[proc] = ttbar_had_xs_model.get_coeff('taggedMass_1btag_000-040tau', proc).get_value(uncertainty_values)
		if proc == "TTBAR":
			print test[proc]
# 		# skip signal processes we are not interested in:
# 		# if proc in ttbar_had_xs_model.signal_processes and proc not in signal_processes[i]: continue
# 		# test[proc] = ttbar_had_xs_model.get_coeff('taggedMass_all_pt', proc).get_value(uncertainty_values)
# 	# print parVals
# # 
	# histos = evaluate_prediction(ttbar_had_xs_model, parameter_values, include_signal = False)
	# # histos = evaluate_prediction(ttbar_had_xs_model, parameter_values, include_signal = True) ## ttbar signal
	# write_histograms_to_rootfile(histos, 'histos-'+filename)

	## Confidence intervals
	res = pl_interval(ttbar_had_xs_model, input = 'data', n=1,cls = [cl_1sigma, cl_2sigma], signal_process_groups = signal_process_groups,parameter = 'ttbar_rate')
	print res


	# model_summary(ttbar_had_xs_model,lnmode='1sigma')
	# report.write_html('htmlout')




