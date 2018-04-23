# Check which sequences in HiMAP db have multiple phyla
# and remove them (usually very small number of strains.)
#

library(himap)
# Load all from himap
library(data.table)

# Load V3-V4 reference copy table
ref_filename = system.file(
      'database',
      himap_option('blast_dbs')[Hypervariable_region=='V4', table],
      package='himap'
   )

ref.dt = fread(ref_filename, colClasses=c('character', 'character', 'integer'))
ref.dt[, genus := sub('^([^_]+)_.*', '\\1', strain_name)]
ref.dt[, strain_name_norrn := sub('_@rrn[0-9]+', '', strain_name)]

# Load taxonomy
tax.dt = load_taxonomy()
tax_g.dt = unique(tax.dt[, .(phylum, class, order, family, genus)])

# Merge tables
ref_tax.dt = merge(ref.dt, tax_g.dt, by='genus')

# Get all variant_ids with multiple phyla
multi_phyla.dt = ref_tax.dt[, if(length(unique(phylum)) > 1) .SD, by=variant_id]

# Now get only variant_ids for which 99% or more strains have only 1 phylum
multi_phyla_weird.dt = multi_phyla.dt[, {
   ft = sort(table(phylum), decreasing=T)
   list('major_phylum'=names(ft)[1], 'coverage_n'=unname(ft[1]),
        'total'=sum(ft),
        'coverage_pct'=100*unname(ft[1])/sum(ft))
}, by=variant_id]

variant_ids_fix = multi_phyla_weird.dt[coverage_pct>=99, variant_id]
major_phyla_fix = multi_phyla_weird.dt[coverage_pct>=99, major_phylum]
strain_names = sub('_@rrn[0-9]+', '', unname(unlist(mapply(function (phy, vid) {
   ref_tax.dt[phylum != phy & variant_id == vid, strain_name]
}, major_phyla_fix, variant_ids_fix))))

ref_fix.dt = ref.dt[!(strain_name_norrn %in% strain_names), 1:3]

write_table(ref_fix.dt, ref_filename)
