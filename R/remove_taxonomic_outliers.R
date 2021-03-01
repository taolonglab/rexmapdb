# Check which sequences in HiMAP db have multiple ranks
# and remove them (usually very small number of strains.)
# if they are < 1% in that rank.

library(rexmap)
library(data.table)
library(optparse)
library(parallel)

# Input and output file (gets overwritten)
# args = commandArgs(trailingOnly=TRUE)
# ref_filename = args[1]
# out_filename = args[2]
# out_exclusions = args[3]
# cov_pct_th = 99

# DEBUG
# setwd('~/cloud/research/microbiome/rexmapdb/')
# ref_filename = 'data/V1-V2_V3-V4_V5-V6_V7-V8_hang22_2021-02-13_wrefseq_table_unique_variants.txt'
# out_filename = 'data/V1-V2_V3-V4_V5-V6_V7-V8_hang22_2021-02-13_wrefseq_table_unique_variants_R.txt'
# ref_fasta = 'data/V1-V2_V3-V4_V5-V6_V7-V8_hang22_2021-02-13_wrefseq_sequences_unique_variants.fasta'
# out_fasta = 'data/V1-V2_V3-V4_V5-V6_V7-V8_hang22_2021-02-13_wrefseq_sequences_unique_variants_R.fasta'
# cov_pct_th = 99.0
# verbose = TRUE

option_list = list(
  make_option(c('--input-table'), type='character',
              help='Input wrefseq_table_unique_variants table.'),
  make_option(c('--input-fasta'), type='character',
              help='Input wrefseq_sequences_unique_variants.fasta.'),
  make_option(c('--output-table'), type='character',
              help='Output filtered table.'),
  make_option(c('--output-fasta'), type='character',
              help='Output filtered FASTA.'),
  make_option(c('--output-excl'), type='character',
              help='Output list of filtered excluded strains.'),
  make_option(c('--cov-pct-th'), type='double', default=99.0,
              help='Coverage % threshold. Exclude sequence variants from strains
               that cover this % of specific taxonomic rank. E.g. if this is set
              to 99, and one sequence variant maps to 99 Bacteroides strains and
              1 E. coli strain, E. coli strain will be excluded.'),
  make_option(c('--verbose'), type='logical', default=TRUE,
              help='Verbose print output.')
)

opt = parse_args(OptionParser(option_list=option_list))

ref_filename = opt$`input-table`
ref_fasta = opt$`input-fasta`
out_filename = opt$`output-table`
out_fasta = opt$`output-fasta`
out_exclusions = opt$`output-excl`
cov_pct_th = opt$`cov-pct-th`
verbose = opt$`verbose`

p = function (...) {
  if (verbose) {
    cat(...)
  }
}

#-------------------------- Load input files -----------------------------------
# Input FASTA
p('* Loading input FASTA...')
ref_fasta = rexmap:::fasta_reader(ref_fasta)
p(' OK.\n')

# Generate a list with sequences as keys and values being a list with
# list(variant_id-copy_number, c(strains, ...))

# For each meta-seq pair with an index i do this
ref_fasta_nreads = length(ref_fasta$meta)
p('* Extract strain names from meta-data.')
ref_fasta_variant_ids = sapply(ref_fasta$meta, function (m)
  sub('^([0-9]+)\\-.*$', '\\1', m), USE.NAMES=F)
p('.')
ref_fasta_strainlabels = sapply(ref_fasta$meta, function (m)
  sub('^([0-9]+)\\-', '', m), USE.NAMES=F)
p('.')
ref_fasta_strainlabels_list = lapply(ref_fasta_strainlabels, function (s)
  strsplit(s, ';', fixed=T))
p(' OK.\n')

# Input table
p('* Loading input table...')
ref.dt = fread(ref_filename, colClasses=c('character', 'character', 'integer',
                                          'character'), showProgress=F)
p(' OK.\n')
p('* Processing table...')
if (names(ref.dt)[2] == 'strain_name' & !('variant_name' %in% names(ref.dt))) {
  names(ref.dt)[2] = 'variant_name'
}
ref.dt[, genus := sub('^([^_]+)_.*', '\\1', variant_name)]
ref.dt[genus=='Candidatus', genus := sub('^[^_]+_([^_]+)_.*', '\\1', variant_name)]
ref.dt[, strain_name_norrn := sub('_@rrn[0-9]+', '', variant_name)]
p(' OK.\n')

# Load taxonomy
p('* Load and organize RExMap taxonomy...')
tax.dt = load_taxonomy()
tax_g.dt = unique(tax.dt[, .(phylum, class, order, family, genus)])
tax_g.dt = rbindlist(list(
  tax_g.dt,
  tax_g.dt[genus=='Rheinheimera', .(genus='Pararheinheimera'),
           by=.(phylum, class, order, family)]
))
unique_genera = tax_g.dt[, unique(genus)]
p(' OK.\n')





p('* Filtering strains:\n')
# Strain names to exclude
p('  - Double strain names...')
strain_names = c()

# Find strain names that contain two species names
strain_names_double = ref.dt[
  strain_name_norrn %like% '^[^_]+_[^_]+_[^_]+_.*' &
  sub('[^_]+_[^_]+_([^_]+)_.*', '\\1', strain_name_norrn) %in% unique_genera,
  strain_name_norrn
]
# Update strain names
strain_names = union(strain_names, strain_names_double)
p(' OK.\n')

p('  - Phylum outliers (', cov_pct_th, '%)...', sep='')
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

variant_ids_fix = multi_phyla_weird.dt[coverage_pct>=cov_pct_th, variant_id]
major_phyla_fix = multi_phyla_weird.dt[coverage_pct>=cov_pct_th, major_phylum]
strain_names_phylum = sub('_@rrn[0-9]+', '', unname(unlist(mapply(function (phy, vid) {
   ref_tax.dt[!(strain_name_norrn %in% strain_names) & phylum != phy & variant_id == vid, variant_name]
}, major_phyla_fix, variant_ids_fix))))

# Update strain names
strain_names = union(strain_names, strain_names_phylum)
p(' OK.\n')


p('  - Class outliers (', cov_pct_th, '%)...', sep='')
# Now look for multiple classes (without multiple phyla strain names)
multi_class.dt = ref_tax.dt[!(strain_name_norrn %in% strain_names),
                            if(length(unique(class)) > 1) .SD, by=variant_id]
multi_class_weird.dt = multi_class.dt[, {
   ft = sort(table(class), decreasing=T)
   list('major_class'=names(ft)[1], 'coverage_n'=unname(ft[1]),
        'total'=sum(ft),
        'coverage_pct'=100*unname(ft[1])/sum(ft))
}, by=variant_id]
variant_ids_fix_class = multi_class_weird.dt[coverage_pct>=cov_pct_th, variant_id]
major_class_fix = multi_class_weird.dt[coverage_pct>=cov_pct_th, major_class]
strain_names_class = sub('_@rrn[0-9]+', '', unname(unlist(mapply(function (cls, vid) {
   ref_tax.dt[!(strain_name_norrn %in% strain_names) & class != cls &
                variant_id == vid, variant_name]
}, major_class_fix, variant_ids_fix_class))))

# Update strain names
strain_names = union(strain_names, strain_names_class)
p(' OK.\n')

p('  - Order outliers (', cov_pct_th, '%)...', sep='')
# Now check order
multi_order.dt = ref_tax.dt[!(strain_name_norrn %in% strain_names),
                            if(length(unique(order)) > 1) .SD, by=variant_id]
multi_order_weird.dt = multi_order.dt[, {
   ft = sort(table(order), decreasing=T)
   list('major_order'=names(ft)[1], 'coverage_n'=unname(ft[1]),
        'total'=sum(ft),
        'coverage_pct'=100*unname(ft[1])/sum(ft))
}, by=variant_id]
variant_ids_fix_order = multi_order_weird.dt[coverage_pct>=cov_pct_th, variant_id]
major_order_fix = multi_order_weird.dt[coverage_pct>=cov_pct_th, major_order]
strain_names_order = sub('_@rrn[0-9]+', '', unname(unlist(mapply(function (cls, vid) {
   ref_tax.dt[!(strain_name_norrn %in% strain_names) & order != cls &
                variant_id == vid, variant_name]
}, major_order_fix, variant_ids_fix_order))))

# Update strain names
strain_names = union(strain_names, strain_names_order)
p(' OK.\n')

p('  - Family outliers (', cov_pct_th, '%)...', sep='')
# Now check family
multi_family.dt = ref_tax.dt[!(strain_name_norrn %in% strain_names),
                            if(length(unique(family)) > 1) .SD, by=variant_id]
multi_family_weird.dt = multi_family.dt[, {
   ft = sort(table(family), decreasing=T)
   list('major_family'=names(ft)[1], 'coverage_n'=unname(ft[1]),
        'total'=sum(ft),
        'coverage_pct'=100*unname(ft[1])/sum(ft))
}, by=variant_id]
variant_ids_fix_family = multi_family_weird.dt[coverage_pct>=cov_pct_th, variant_id]
major_family_fix = multi_family_weird.dt[coverage_pct>=cov_pct_th, major_family]
strain_names_family = sub('_@rrn[0-9]+', '', unname(unlist(mapply(function (cls, vid) {
   ref_tax.dt[!(strain_name_norrn %in% strain_names) & family != cls &
                variant_id == vid, variant_name]
}, major_family_fix, variant_ids_fix_family))))

# Update strain names
strain_names = union(strain_names, strain_names_family)
p(' OK.\n')


p('* Generating exclusion table...')
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
p(' OK.\n')

#------------------- WRITE DATA -----------------------------------
p('* Writing exclusion table...')
# Write table with excluded strains
write.table(excl.dt, out_exclusions, sep='\t', quote=F, row.names=F)
p(' OK.\n')

p('* Writing output table...')
# Write fixed reference table with copy numbers
ref_fix.dt = ref.dt[!(strain_name_norrn %in% strain_names), 1:4]
write_table(ref_fix.dt[, 1:3], out_filename, sep='\t')
p(' OK.\n')

# Write filtered FASTA file
p('* Generating new FASTA meta-data labels...')
ref_fasta_new.dt = ref_fix.dt[
  , .(strain_names=paste(unique(sort(sub('_@rrn[0-9]+', '', variant_name))),
                         collapse=','))
  , by=.(variant_id, copy_number, seq)][
    , .(meta=paste0(.BY[[1]], '-',
             paste(copy_number, paste('(', strain_names, ')', sep=''),
                   sep=':',
                   collapse=';')
             ))
    , by=.(variant_id, seq)]
p(' OK.\n')

p('* Writing output FASTA...')
rexmap:::fasta_writer(meta=ref_fasta_new.dt[, meta],
                      seqs=ref_fasta_new.dt[, seq],
                      output=out_fasta)
p(' OK.\n')


