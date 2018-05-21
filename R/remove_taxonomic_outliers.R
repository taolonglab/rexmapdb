# Check which sequences in HiMAP db have multiple ranks
# and remove them (usually very small number of strains.)
# if they are < 1% in that rank.

library(himap)
library(data.table)

# Input and output file (gets overwritten)
args = commandArgs(trailingOnly=TRUE)
ref_filename = args[1]
out_filename = args[2]
out_exclusions = args[3]

# Load V3-V4 reference copy table
# ref_filename = system.file(
#       'database',
#       himap_option('blast_dbs')[Hypervariable_region=='V4', table],
#       package='himap'
#    )

# ref_filename = '~/cloud/research/microbiome/genomes/data/vregions_db/V3-V4_337F-805R_hang21_wrefseq_table_unique_variants.txt'
# ref_filename = '~/cloud/research/microbiome/genomes/data/vregions_db/V4_515F-805R_hang22_wrefseq_table_unique_variants.txt'
#
# out_filename = '~/cloud/research/microbiome/himap/inst/database/V3-V4_337F-805R_hang21_wrefseq_table_unique_variants_R.txt'
# out_filename = '~/cloud/research/microbiome/himap/inst/database/V4_515F-805R_hang22_wrefseq_table_unique_variants_R.txt'


ref.dt = fread(ref_filename, colClasses=c('character', 'character', 'integer',
                                          'character'))
ref.dt[, genus := sub('^([^_]+)_.*', '\\1', strain_name)]
ref.dt[genus=='Candidatus', genus := sub('^[^_]+_([^_]+)_.*', '\\1', strain_name)]
ref.dt[, strain_name_norrn := sub('_@rrn[0-9]+', '', strain_name)]

# Load taxonomy
tax.dt = load_taxonomy()
tax_g.dt = unique(tax.dt[, .(phylum, class, order, family, genus)])
tax_g.dt = rbindlist(list(
  tax_g.dt,
  tax_g.dt[genus=='Rheinheimera', .(genus='Pararheinheimera'),
           by=.(phylum, class, order, family)]
))
unique_genera = tax_g.dt[, unique(genus)]

# Strain names to exclude
strain_names = c()

# Find strain names that contain two species names
strain_names_double = ref.dt[
  strain_name_norrn %like% '^[^_]+_[^_]+_[^_]+_.*' &
  sub('[^_]+_[^_]+_([^_]+)_.*', '\\1', strain_name_norrn) %in% unique_genera,
  strain_name_norrn
]

# Update strain names
strain_names = union(strain_names, strain_names_double)

# Merge tables
ref_tax.dt = merge(ref.dt, tax_g.dt, by='genus')

# Get all variant_ids with multiple phyla
multi_phyla.dt = ref_tax.dt[!(strain_name_norrn %in% strain_names),
                            if(length(unique(phylum)) > 1) .SD, by=variant_id]

# Now get only variant_ids for which 99% or more strains have only 1 phylum
multi_phyla_weird.dt = multi_phyla.dt[, {
   ft = sort(table(phylum), decreasing=T)
   list('major_phylum'=names(ft)[1], 'coverage_n'=unname(ft[1]),
        'total'=sum(ft),
        'coverage_pct'=100*unname(ft[1])/sum(ft))
}, by=variant_id]

variant_ids_fix = multi_phyla_weird.dt[coverage_pct>=99, variant_id]
major_phyla_fix = multi_phyla_weird.dt[coverage_pct>=99, major_phylum]
strain_names_phylum = sub('_@rrn[0-9]+', '', unname(unlist(mapply(function (phy, vid) {
   ref_tax.dt[!(strain_name_norrn %in% strain_names) & phylum != phy & variant_id == vid, strain_name]
}, major_phyla_fix, variant_ids_fix))))

# Update strain names
strain_names = union(strain_names, strain_names_phylum)

# Now look for multiple classes (without multiple phyla strain names)
multi_class.dt = ref_tax.dt[!(strain_name_norrn %in% strain_names),
                            if(length(unique(class)) > 1) .SD, by=variant_id]
multi_class_weird.dt = multi_class.dt[, {
   ft = sort(table(class), decreasing=T)
   list('major_class'=names(ft)[1], 'coverage_n'=unname(ft[1]),
        'total'=sum(ft),
        'coverage_pct'=100*unname(ft[1])/sum(ft))
}, by=variant_id]
variant_ids_fix_class = multi_class_weird.dt[coverage_pct>=99, variant_id]
major_class_fix = multi_class_weird.dt[coverage_pct>=99, major_class]
strain_names_class = sub('_@rrn[0-9]+', '', unname(unlist(mapply(function (cls, vid) {
   ref_tax.dt[!(strain_name_norrn %in% strain_names) & class != cls &
                variant_id == vid, strain_name]
}, major_class_fix, variant_ids_fix_class))))

# Update strain names
strain_names = union(strain_names, strain_names_class)

# Now check order
multi_order.dt = ref_tax.dt[!(strain_name_norrn %in% strain_names),
                            if(length(unique(order)) > 1) .SD, by=variant_id]
multi_order_weird.dt = multi_order.dt[, {
   ft = sort(table(order), decreasing=T)
   list('major_order'=names(ft)[1], 'coverage_n'=unname(ft[1]),
        'total'=sum(ft),
        'coverage_pct'=100*unname(ft[1])/sum(ft))
}, by=variant_id]
variant_ids_fix_order = multi_order_weird.dt[coverage_pct>=99, variant_id]
major_order_fix = multi_order_weird.dt[coverage_pct>=99, major_order]
strain_names_order = sub('_@rrn[0-9]+', '', unname(unlist(mapply(function (cls, vid) {
   ref_tax.dt[!(strain_name_norrn %in% strain_names) & order != cls &
                variant_id == vid, strain_name]
}, major_order_fix, variant_ids_fix_order))))

# Update strain names
strain_names = union(strain_names, strain_names_order)

# Now check family
multi_family.dt = ref_tax.dt[!(strain_name_norrn %in% strain_names),
                            if(length(unique(family)) > 1) .SD, by=variant_id]
multi_family_weird.dt = multi_family.dt[, {
   ft = sort(table(family), decreasing=T)
   list('major_family'=names(ft)[1], 'coverage_n'=unname(ft[1]),
        'total'=sum(ft),
        'coverage_pct'=100*unname(ft[1])/sum(ft))
}, by=variant_id]
variant_ids_fix_family = multi_family_weird.dt[coverage_pct>=99, variant_id]
major_family_fix = multi_family_weird.dt[coverage_pct>=99, major_family]
strain_names_family = sub('_@rrn[0-9]+', '', unname(unlist(mapply(function (cls, vid) {
   ref_tax.dt[!(strain_name_norrn %in% strain_names) & family != cls &
                variant_id == vid, strain_name]
}, major_family_fix, variant_ids_fix_family))))

# Update strain names
strain_names = union(strain_names, strain_names_family)

# Generate a table with exclusions, including the reason
excl.dt = unique(data.table(
  strain = c(strain_names_double, strain_names_phylum, strain_names_class,
             strain_names_order, strain_names_family),
  reason = factor(
    c(rep('double species', length(strain_names_double)),
      rep('phylum outlier', length(strain_names_phylum)),
      rep('class outlier', length(strain_names_class)),
      rep('order outlier', length(strain_names_order)),
      rep('family outlier', length(strain_names_family))
    ),
    levels=c('double species', 'phylum outlier', 'class outlier',
             'order outlier', 'family outlier')
  )
))
setorder(excl.dt, reason, strain)

#------------------- WRITE DATA -----------------------------------
# Write table with excluded strains
write.table(excl.dt, out_exclusions, sep='\t', quote=F, row.names=F)

# Write fixed reference table with copy numbers
ref_fix.dt = ref.dt[!(strain_name_norrn %in% strain_names), 1:3]
write_table(ref_fix.dt, out_filename, sep='\t', quote=F, row.names=F)
