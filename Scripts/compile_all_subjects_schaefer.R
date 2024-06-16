output_path = '/project/bbl_roalf_longglucest/sandbox/ally/test_schaefer/output_measures'

for (str in c('INV2', 'UNI')) {
    output = file.path(output_path, str)

        all_data = dplyr::tibble()
        all_inv2 = dplyr::tibble()

        for (sub_dir in list.dirs(output, recursive = FALSE)) {
            subject = basename(sub_dir)
            sub_data = file.path(sub_dir, paste0(subject, 
					     paste0('-Schaefer2018ROI-GluCEST-100P-17N-measures_', str, '.tsv')))
            sub_data = readr::read_tsv(sub_data)
            all_data = dplyr::bind_rows(all_data, sub_data)
        
	    sub_data = file.path(sub_dir, paste0(subject,
	        				 paste0('-Schaefer2018ROI-INV2-100P-17N-measures_', str, '.tsv')))
	    sub_data = readr::read_tsv(sub_data)
	    all_inv2 = dplyr::bind_rows(all_inv2, sub_data)
	
        }

        for (network in c("Cont", "Default", "DorsAttn", "Limbic", "SalVentAttn", "SomMot", "Vis")) {
            all_subs = dplyr::tibble()

            for (sub_dir in list.dirs(output, recursive = FALSE)) {
                subject = basename(sub_dir)
                sub_data = file.path(sub_dir, paste0(subject,
                            paste0('-2d-GluCEST-s100_7-', network, '-measures_', str, '.tsv')))
                sub_data = readr::read_tsv(sub_data)
                all_subs = dplyr::bind_rows(all_subs, sub_data)
            }

            readr::write_tsv(all_subs, file.path(output,
                                paste0('GluCEST-s100_7-', network, '-Measures_', str, '.tsv')))
        }


    readr::write_tsv(all_data, file.path(output, 
				         paste0('GluCEST-Schaefer2018ROI-100P-17N-Measures_', str, '.tsv')))
    readr::write_tsv(all_inv2, file.path(output,
					 paste0(str, '-Schaefer2018ROI-100P-17N-Measures_', str, '.tsv')))
}
