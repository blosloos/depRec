# include <Rcpp.h>

using namespace Rcpp;


// set RT windows for a cluster ...............................................//
std::vector<double> set_RT_window(
	NumericMatrix peaks,
	std::vector<int>& for_peaks,
	double dret
){
	
	std::vector<double> RT_lim(2);
	int size = for_peaks.size(), n;
	double RT_span;
	
	RT_lim[0] = peaks(for_peaks[0] - 1, 2);
	RT_lim[1]  = peaks(for_peaks[0] - 1, 2);
	if(size > 1){
		for(n = 1; n < size; n++){
			if( peaks(for_peaks[n] - 1, 2) < RT_lim[0] ){RT_lim[0] = peaks(for_peaks[n] - 1, 2);};
			if( peaks(for_peaks[n] - 1, 2) > RT_lim[1] ){RT_lim[1] = peaks(for_peaks[n] - 1, 2);};		
		
		}
	}
	RT_span = ((dret - (RT_lim[1] - RT_lim[0])) / 2);
	RT_lim[0] = (RT_lim[0] - RT_span);
	RT_lim[1] = (RT_lim[1] + RT_span);
	
	return RT_lim;
}

// set tolerance windows for a mass ..........................................//
std::vector<double> set_mass_window(
	NumericMatrix peaks,
	std::vector<int>& for_peaks,
	bool ppm,
	double dmass
){
	
	std::vector<double> mass_lim(2);
	int size = for_peaks.size(), n;
	double mass, delmass;
	
	mass = peaks(for_peaks[0] - 1, 0);
	if(ppm){
		delmass = (dmass * mass / 1000000);
	}else{
		delmass = (dmass / 1000);
	};
	mass_lim[0] = mass - delmass;
	mass_lim[1] = mass + delmass;	
	if(size > 1){
		for(n = 1; n < size; n++){
			mass = peaks(for_peaks[n] - 1, 0);
			if(ppm){
				delmass = (dmass * mass / 1000000);
			}else{
				delmass = (dmass / 1000);
			};	
			if( (mass + delmass) < mass_lim[1] ){mass_lim[1] = (mass + delmass);};
			if( (mass - delmass) > mass_lim[0] ){mass_lim[0] = (mass - delmass);};				
			
		}
	}
	
	return mass_lim;
}

bool check_fits_cluster(
	NumericMatrix peaks,
	std::vector<int>& for_peaks,	
	std::vector<double>& mass_lim_low,
	std::vector<double>& mass_lim_up,	
	std::vector<double>& RT_lim_low,
	std::vector<double>& RT_lim_up,	
	int which_cluster,
	std::vector<int>& x_new,
	bool ppm,
	double dmass	
){
	int m, new_size = x_new.size();
	
	// new_peak ranges within RT limits?
	for(m = 0; m < new_size; m++){
		if(peaks(x_new[m] - 1, 2) < RT_lim_low[which_cluster]) return false;
		if(peaks(x_new[m] - 1, 2) > RT_lim_up[which_cluster]) return false;	
	}
		
	// new_peak overlaps with mass limits?	
	double mass, delmass;	
	for(m = 0; m < new_size; m++){
		mass = peaks(x_new[m] - 1, 0);
		if(ppm){
			delmass = (dmass * mass / 1000000);
		}else{
			delmass = (dmass / 1000);
		};
		if((mass + delmass) < mass_lim_low[which_cluster]) return false;
		if((mass - delmass) > mass_lim_up[which_cluster]) return false;
	}
	
	// already a peak with same sampleID in the cluster?
	int n, size = for_peaks.size(), sample_ID;
	for(m = 0; m < new_size; m++){	
		sample_ID = peaks(x_new[m] - 1, 3);
		for(n = 0; n < size; n++){	
			if(peaks(for_peaks[n] - 1, 3) == sample_ID) return false;
		}
	}
	
	return true;
}

// [[Rcpp::export]]
IntegerVector extractProfiles(
	NumericMatrix peaks,
	IntegerVector in_order,
	double dmass,
	bool ppm,
	double dret
){

	unsigned int n_rows = peaks.nrow(), n = 0, m, at_peak, which_cluster;
	bool fits_to_cluster;
	std::vector<int> x_new(1);
	std::vector<int> x_use;
	std::vector<double> mass_lim_low;
	std::vector<double> mass_lim_up;
	std::vector<double> RT_lim_low;
	std::vector<double> RT_lim_up;
	std::vector<double> some_values(2);	
	std::vector< std::vector<int> > cluster_with_peaks; 
	
	
	// initiate first cluster .................................................... //
	at_peak = in_order[n];
	x_new[0] = at_peak;
	cluster_with_peaks.push_back(x_new); 
	some_values = set_RT_window(peaks, x_new, dret);
	RT_lim_low.push_back(some_values[0]);
	RT_lim_up.push_back(some_values[1]);
	some_values = set_mass_window(peaks, x_new, ppm, dmass);	
	mass_lim_low.push_back(some_values[0]);
	mass_lim_up.push_back(some_values[1]);	
	// cluster over remaining peaks .............................................. //
	if(n_rows > 1){
		for(n = 1; n < n_rows; n++){
			at_peak = in_order[n];		
			x_new[0] = at_peak;
			fits_to_cluster = false; // check: fits to any existing cluster? ......... //
			for(m = 0; m < cluster_with_peaks.size(); m++){
				if( check_fits_cluster(peaks, cluster_with_peaks[m], mass_lim_low, mass_lim_up,	RT_lim_low, RT_lim_up, m, x_new, ppm, dmass) ) {
					fits_to_cluster = true;
					which_cluster = m;
					break;
				}
			}
			if(fits_to_cluster){ // add to existing cluster 
				x_use = cluster_with_peaks[which_cluster];
				x_use.push_back(at_peak);
				cluster_with_peaks[which_cluster] = x_use;
				some_values = set_RT_window(peaks, x_use, dret);
				RT_lim_low[which_cluster] = some_values[0];
				RT_lim_up[which_cluster] = some_values[1];				
				some_values = set_mass_window(peaks, x_use, ppm, dmass);	
				mass_lim_low[which_cluster] = some_values[0];
				mass_lim_up[which_cluster] = some_values[1];				
			}else{ // create new cluster		
				cluster_with_peaks.push_back(x_new); 
				some_values = set_RT_window(peaks, x_new, dret);
				RT_lim_low.push_back(some_values[0]);
				RT_lim_up.push_back(some_values[1]);
				some_values = set_mass_window(peaks, x_new, ppm, dmass);	
				mass_lim_low.push_back(some_values[0]);
				mass_lim_up.push_back(some_values[1]);				
			}
		}
	}
	
	// get some details back ... use "List" as return value instead!
	//List results(5); // [0]: clusters_with_peaks; [1]: mass_lim_low; [2]: mass_lim_up; [3]: RT_lim_low; [4]: RT_lim_up;
	//results[0] = cluster_with_peaks;
	//results[1] = mass_lim_low;
	//results[2] = mass_lim_up;	
	//results[3] = RT_lim_low;
	//results[4] = RT_lim_up;
	//return results;
	
	IntegerVector peaks_in_cluster(n_rows);
	for(m = 0; m < cluster_with_peaks.size(); m++){
		x_use = cluster_with_peaks[m];
		for(n = 0; n < x_use.size(); n++){
			peaks_in_cluster[(x_use[n] - 1)] = (m + 1);
		}
	}
	
	return peaks_in_cluster;
	
}

// [[Rcpp::export]]
IntegerVector extractProfiles_replicates(
	NumericMatrix peaks,
	IntegerVector in_order,
	double dmass,
	bool ppm,
	double dret,
	IntegerVector pregroup
){

	unsigned int n_rows = peaks.nrow(), n = 0, m, at_peak, which_cluster;
	int replic_group;
	bool fits_to_cluster;
	std::vector<int> x_new(1);
	std::vector<int> x_use;
	std::vector<bool> done_group;
	done_group.assign (n_rows, false);
	std::vector<double> mass_lim_low;
	std::vector<double> mass_lim_up;
	std::vector<double> RT_lim_low;
	std::vector<double> RT_lim_up;
	std::vector<double> some_values(2);	
	std::vector< std::vector<int> > cluster_with_peaks; 
	
	
	// initiate first cluster .................................................... //
	at_peak = in_order[n];
	done_group[at_peak - 1] = true;
	x_new[0] = at_peak;
	replic_group = pregroup[at_peak - 1];
	if(replic_group != 0){
		for(m = 0; m < n_rows; m++){
			if((pregroup[m] == replic_group) & (m != (at_peak -1))){
				done_group[m] = true;
				x_new.push_back(m + 1);
			}
		}
	}
	cluster_with_peaks.push_back(x_new); 
	some_values = set_RT_window(peaks, x_new, dret);
	RT_lim_low.push_back(some_values[0]);
	RT_lim_up.push_back(some_values[1]);
	some_values = set_mass_window(peaks, x_new, ppm, dmass);	
	mass_lim_low.push_back(some_values[0]);
	mass_lim_up.push_back(some_values[1]);	
	// cluster over remaining peaks .............................................. //
	if(n_rows > 1){
		for(n = 1; n < n_rows; n++){
			at_peak = in_order[n];	
			if(done_group[at_peak - 1]) continue; // assigned as part of a replicate group?			
			x_new.resize(1);
			x_new[0] = at_peak;			
			replic_group = pregroup[at_peak - 1];
			if(replic_group != 0){
				for(m = 0; m < n_rows; m++){
					if((pregroup[m] == replic_group) & (m != (at_peak -1))){
						done_group[m] = true;
						x_new.push_back(m + 1);
					}
				}
			}
			fits_to_cluster = false; // check: fits to any existing cluster? ......... //
			for(m = 0; m < cluster_with_peaks.size(); m++){
				if( check_fits_cluster(peaks, cluster_with_peaks[m], mass_lim_low, mass_lim_up,	RT_lim_low, RT_lim_up, m, x_new, ppm, dmass) ) {
					fits_to_cluster = true;
					which_cluster = m;
					break;
				}
			}
			if(fits_to_cluster){ // add to existing cluster 
				x_use = cluster_with_peaks[which_cluster];
				x_use.reserve(x_use.size() + x_new.size());
				x_use.insert(x_use.end(), x_new.begin(), x_new.end());
				cluster_with_peaks[which_cluster] = x_use;
				some_values = set_RT_window(peaks, x_use, dret);
				RT_lim_low[which_cluster] = some_values[0];
				RT_lim_up[which_cluster] = some_values[1];				
				some_values = set_mass_window(peaks, x_use, ppm, dmass);	
				mass_lim_low[which_cluster] = some_values[0];
				mass_lim_up[which_cluster] = some_values[1];				
			}else{ // create new cluster		
				cluster_with_peaks.push_back(x_new); 
				some_values = set_RT_window(peaks, x_new, dret);
				RT_lim_low.push_back(some_values[0]);
				RT_lim_up.push_back(some_values[1]);
				some_values = set_mass_window(peaks, x_new, ppm, dmass);	
				mass_lim_low.push_back(some_values[0]);
				mass_lim_up.push_back(some_values[1]);				
			}
		}
	}
	
	// get some details back ... use "List" as return value instead!
	//List results(5); // [0]: clusters_with_peaks; [1]: mass_lim_low; [2]: mass_lim_up; [3]: RT_lim_low; [4]: RT_lim_up;
	//results[0] = cluster_with_peaks;
	//results[1] = mass_lim_low;
	//results[2] = mass_lim_up;	
	//results[3] = RT_lim_low;
	//results[4] = RT_lim_up;
	//return results;
	
	IntegerVector peaks_in_cluster(n_rows);
	for(m = 0; m < cluster_with_peaks.size(); m++){
		x_use = cluster_with_peaks[m];
		for(n = 0; n < x_use.size(); n++){
			peaks_in_cluster[(x_use[n] - 1)] = (m + 1);
		}
	}
	
	return peaks_in_cluster;
	
}







