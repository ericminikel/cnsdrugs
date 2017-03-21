
options(stringsAsFactors=FALSE)

# from http://zinc.docking.org/catalogs/fda
zinc_fda_prop = read.table('c:/sci/037cinf/data/EPA-DSSTOX/fda_prop.xls',header=TRUE,sep='\t',quote="",comment.char="")
# unzipped from http://www.epa.gov/ncct/dsstox/StructureDataFiles/FDAMDD_DownloadFiles/FDAMDD_v3b_1216_15Feb2008.zip and exported to plain text
fda_annotations = read.table('c:/sci/037cinf/data/EPA-DSSTOX/FDAMDD_v3b_1216_15Feb2008_nostructures.txt',header=TRUE,sep='\t',quote="",comment.char="")

dim(fda_annotations) # 1216 34
dim(zinc_fda_prop) # 3180 14

# can we join on SMILES?
sum(fda_annotations$STRUCTURE_SMILES %in% zinc_fda_prop$SMILES) # 2
sum(fda_annotations$STRUCTURE_Parent_SMILES %in% zinc_fda_prop$SMILES) # 4
# not a case issue:
sum(toupper(fda_annotations$STRUCTURE_Parent_SMILES) %in% toupper(zinc_fda_prop$SMILES)) # 4

# can we join on molecular weight?
sum(fda_annotations$STRUCTURE_MolecularWeight %in% zinc_fda_prop$MWT) #2
# ZINC rounds to 3 decimal places, EPA to only 4. 
sum(round(fda_annotations$STRUCTURE_MolecularWeight,3) %in% 
      round(zinc_fda_prop$MWT,3)) # 24 matches when rounded to 3 decimal places
sum(round(fda_annotations$STRUCTURE_MolecularWeight,2) %in% 
      round(zinc_fda_prop$MWT,2)) # 302 matches when rounded to 2 decimal places
sum(round(fda_annotations$STRUCTURE_MolecularWeight,1) %in% 
      round(zinc_fda_prop$MWT,1)) # 768 matches when rounded to 1 decimal place

# but rounding to 1 decimal place means the match is no longer unique b/c 860 < 1216
length(unique(round(fda_annotations$STRUCTURE_MolecularWeight,1))) # 860




names(fda_annotations)

length(unique(fda_annotations$TherapeuticCategory)) # 319

# https://www.unitedhealthcareonline.com/ccmcontent/ProviderII/UHC/en-US/Assets/ProviderStaticFiles/ProviderStaticFilesPdf/Tools%20and%20Resources/Pharmacy%20Resources/PDL_Phys_Bk.pdf
uhc_cns_drug_list = read.table('c:/sci/037cinf/data/uhc-cns-drugs-pdf-list.txt',header=FALSE,sep="|") # read file, choose a delimiter absent from the file
cns_drug_string = paste(uhc_cns_drug_list$V1, sep=" ", collapse=" ") # collapse it into a single huge string
cns_drug_string = tolower(cns_drug_string) # convert to lowercase

fda_annotations$cns_drug = FALSE # initialize a column for whether a cns drug or not
for (i in 1:dim(fda_annotations)[1]) { # loop through b/c grep doesn't support vector operations.
  fda_annotations$cns_drug[i] = length(grep(tolower(fda_annotations$TestSubstance_ChemicalName[i]),cns_drug_string)) > 0 
}
which(fda_annotations$cns_drug) # list of indices
fda_annotations$TestSubstance_ChemicalName[fda_annotations$cns_drug] # list of names

# 122 hits but aripiprazole (abilify) is conspicuously absent.

# NDC is useless:
ndc_product = read.table('c:/sci/037cinf/data/ndc/product.txt',sep='\t',header=TRUE,quote="",comment.char="")
dim(ndc_product) # 71344 18
names(ndc_product)
length(unique(ndc_product$NONPROPRIETARYNAME)) # 12672
length(unique(ndc_product$SUBSTANCENAME)) # 6657
head(ndc_product$SUBSTANCENAME)

# in read.table, choose a delimiter absent from the file such as | in order to just read all as one column.
uhc_cns_drug_list      = read.table('c:/sci/037cinf/data/uhc-cns-drugs-pdf-list.txt',header=FALSE,sep="|") # # https://www.unitedhealthcareonline.com/ccmcontent/ProviderII/UHC/en-US/Assets/ProviderStaticFiles/ProviderStaticFilesPdf/Tools%20and%20Resources/Pharmacy%20Resources/PDL_Phys_Bk.pdf 
humana_cns_drug_list   = read.table('c:/sci/037cinf/data/humana_2014_wi.txt',header=FALSE,sep="|") # http://apps.humana.com/marketing/documents.asp?file=2168465
ha_cns_drug_list       = read.table('c:/sci/037cinf/data/healthalliance-2013.txt',header=FALSE,sep="|") # https://www.healthalliance.org/media/HAMP_formulary.pdf%E2%80%8E
bcbst_cns_drug_list    = read.table('c:/sci/037cinf/data/bcbst-2014.txt',header=FALSE,sep="|") # http://www.bcbstx.com/pdf/drug_list.pdf
caremark_cns_drug_list = read.table('c:/sci/037cinf/data/caremark-oct2013.txt',header=FALSE,sep="|")  # http://www.caremark.com/portal/asset/caremark_drug_list.pdf
all_cns_drug_lists = rbind(uhc_cns_drug_list, humana_cns_drug_list, ha_cns_drug_list, bcbst_cns_drug_list, caremark_cns_drug_list) # concatenate all together
cns_drug_string = paste(all_cns_drug_lists$V1, sep=" ", collapse=" ") # collapse it into a single huge string
cns_drug_string = tolower(cns_drug_string) # convert to lowercase
cns_drug_string

drugbank = read.table('c:/sci/037cinf/data/drugbank.txt/drugbank.5col.pipe',sep='|',header=FALSE,quote="",comment.char="")
colnames(drugbank) = c("generic_name","brand_names","first_group","smiles","indication")

# how to remove non-alphanumeric chars:
# http://stackoverflow.com/questions/8959243/r-remove-non-alphanumeric-symbols-from-a-string

drugbank$cns_drug = FALSE # initialize a column for whether a cns drug or not
for (i in 1:dim(drugbank)[1]) { # loop all drugs.  grep doesn't support vector operations.
  # check for generic names
  generic_name_present = FALSE # initialize
  # some generic names are chemical names with []{}() characters which mess up regexes
  # and obviously aren't called by that name in any insurance plans anyway, so skip those
  if(!length(grep("\\}|\\{|\\[|\\]|\\(|\\)",drugbank$generic_name[i])) > 0) {
    # if grep finds a match between the DrugBank generic name and the cns drug list, set TRUE
    generic_name_present = length(grep(paste("\\b",tolower(drugbank$generic_name[i]),"\\b",sep=""),cns_drug_string)) > 0 
  }
  # check for brand names
  brand_name_present = FALSE # initizalize
  brand_name_vector = strsplit(drugbank$brand_names[i],split=',')[[1]] # split the list of brand names on comma
  # remove empty and "Not Available" entries
  brand_name_vector = brand_name_vector[-which(brand_name_vector=="" | brand_name_vector=="Not Available")]
  # loop through all brand names
  for (brand_name in brand_name_vector) { 
    # some brand names in DrugBank are listed with characters like [] that mess up regexes, so:
    brand_name_alphanumeric = gsub("[^[:alnum:] ]", "", brand_name) # remove nonalphanumeric characters
    # now check if it's in the CNS drug list string
    if (length(grep(paste("\\b",tolower(brand_name_alphanumeric),"\\b",sep=""),cns_drug_string)) > 0) {
      brand_name_present = TRUE
    }
  }
  # either a generic or brand name hit counts for setting cns_drug to TRUE
  drugbank$cns_drug[i] = generic_name_present | brand_name_present
}

which(drugbank$cns_drug) # list of indices
drugbank$generic_name[drugbank$cns_drug & drugbank$first_group == 'Approved'] # list of names

# \\b in the grep pattern is for whole word matching
# http://stackoverflow.com/questions/12885257/r-grep-whole-words-separated-by-special-characters
# if you don't do that, you get:
drugbank[drugbank$generic_name=='Iron',]
# upon looking in the UHC list, this showed up because iron is a substring of buspirone, an anxiolytic

# 155 hits with only the UHC list to grep against.
# 175 with all insurances
# 221 with all insurances and brand name + generic

drugbank_approved = subset(drugbank, first_group=='Approved')

drug_list = drugbank[drugbank$first_group == 'Approved', c('generic_name','cns_drug','smiles')]
drug_list_ordered = drug_list[with(drug_list,order(generic_name)),]
dim(drug_list_ordered) # 1691 3
write.table(drug_list_ordered,'drugs.txt',sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)

cns_list = sort(drugbank$generic_name[drugbank$cns_drug & drugbank$first_group=='Approved'])
write.table(cns_list,'cnsdrugs.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)

