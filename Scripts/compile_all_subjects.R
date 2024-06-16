output_path = '/project/bbl_roalf_longglucest/sandbox/ally/test_get_fdata/output_measures'

for (str in c('INV2', 'UNI')) {
    output = file.path(output_path, str)

    for (type in c("cort", "sub")) {
        all_data = dplyr::tibble()
        all_inv2 = dplyr::tibble()

        for (sub_dir in list.dirs(output, recursive = FALSE)) {
            subject = basename(sub_dir)
            sub_data = file.path(sub_dir, paste0(subject, 
					     paste0('-HarvardOxfordROI-GluCEST-', type, '-measures_', str, '.tsv')))
            sub_data = readr::read_tsv(sub_data)
            all_data = dplyr::bind_rows(all_data, sub_data)
        
	    sub_data = file.path(sub_dir, paste0(subject,
	        				 paste0('-HarvardOxfordROI-INV2-', type, '-measures_', str, '.tsv')))
	    sub_data = readr::read_tsv(sub_data)
	    all_inv2 = dplyr::bind_rows(all_inv2, sub_data)
	
        }
    readr::write_tsv(all_data, file.path(output, 
				         paste0('GluCEST-HarvardOxfordROI-', type, '-Measures_', str, '.tsv')))
    readr::write_tsv(all_inv2, file.path(output,
					 paste0(str, '-HarvardOxfordROI-', type, '-Measures_', str, '.tsv')))
    }
}
