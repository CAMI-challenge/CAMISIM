# -*- coding: utf-8 -*-
import sys, glob, argparse
from TaxonomyNcbi import TaxonomyNcbi


__author__ = 'jessika'

"""
    Step 4: Novelty_category determination (Jessika)

    NAME Novelty.py

    FILE /net/metagenomics/projects/cami_2014/reference_data_preparation/Novelty.py

    DESCRIPTION

    Define the taxonomic 'novelty category' for a list of genomes with given NCBI_IDs (taxon IDs, extracted from the metafile.csv) relative to a list of reference strain sequences (extracted from a directory which includes NCBI reference sequences named [taxID].[nr].fna). The taxon ID reflects the lowest rank until which the genome could be placed in the reference taxonomy, meaning that if the rank is above strain level, it is likely new up to the rank below where the taxon ID belongs to. The 'novelty_category' specifies the highest taxonomic rank that a sequenced strain belongs to for which there is no sequenced genome yet in the list of reference strain sequences. This can be even a higher rank than one minus the rank of the taxon ID, as not all taxa which are known have sequenced genomes available. The output is ADD DeTAILS.


    Algorithm:

    if taxon is a strain
        is included  in the reference?  return 'None': return new_strain

    for each rank  from the taxons rank up to rank 'superkingdom':
         is a taxon from the reference a child of the taxon? return new_(rank -1) // this means that there is a sequenced genomes from this taxon in the reference, and one rank below the rank is new (as apparently not defined in reference taxonomy yet and detectable by us using 16S analysis)

    return new_'superkingdom'

    TODO
        None

    Usage: Novelty.py -ref [directory including ncbi reference sequences] -db [ncbi database] -meta [metafile including column NCBI_ID] -o [output file]

"""

class Novelty():

    def __init__(self, db):
        """
            @param db: usually file named "ncbitax_sqlite.db"
        """
        self.tax = TaxonomyNcbi(db)
        self.ranks = ['strain','species','genus','family','order','class','phylum','superkingdom','root']
        #wrapping:
        self.included_parents_at_rank = dict()
        for rank in self.ranks:
            self.included_parents_at_rank[rank] = set()
        self.included_parents_at_rank['no rank'] = set()   #includes no ranks and all strain ids
        self.included_parents_at_rank['root'].add(1)

    def read_reference(self, ref_dir,excluded=None):
        """
            extracts all taxonomic IDs from the reference directory and stores all contained IDs for each rank
            @param referenceIdsSet:   including reference IDs
        """

        print "extracting included parents at each rank for each reference ID.. This may take a while."

        refernceIdsSet = self.get_taxonomic_ids_from_directory(ref_dir)

        counter = 0
        for ncbiid in refernceIdsSet:
            if excluded is not None and ncbiid in excluded:
                continue

            try:
                ncbiid = int(ncbiid.split("_")[0])
            except:
                ncbiid = ncbiid

            counter += 1
            if counter == 30:
                sys.stdout.write('...')   #tell user script is alive
                sys.stdout.flush()
                counter = 0

            parentSet = self.tax.getParentsNcbidSet(ncbiid)
            parentSet.add(ncbiid)                       #this ncbiid is included too

            for p in parentSet:
                r = self.tax.getRank(p)
                if r not in self.ranks:
                    r = 'no rank'
                elif p is 1:
                    r = 'root'
                self.included_parents_at_rank[r].add(int(p))

        print "reference processing done."

    def compute_novelty_for_metafile(self,metafile,output):
        """
            computes the novelty_category for each NCBI ID in the metafile and updates it to the output file
            (Note that the metafile must include a header with column name 'NCBI_ID'
                                  whereas novelty_category is added if it does not exist)
            @param metafile: usually file named 'metadata_table_[version].csv'#
            @param output:  file for the output
        """

        print "processing information from metafile: " + metafile
        try:
            meta = open(metafile)
        except:
            sys.exit("can not open metafile")
        lines = meta.readlines()
        meta.close()

        header = lines[0]

        id_index = header.split("\t").index("NCBI_ID")
        add_column_for_novelty = False
        try:
            novelty_index = header.split("\t").index("novelty_category")
        #add a column to the metafile if "novelty_category" is not included yet:
        except:
            add_column_for_novelty = True
            header = header.rstrip("\n")
            header += "\tnovelty_category"
            novelty_index = header.split("\t").index("novelty_category")

        try:
            out = open(output,"w")
        except:
            sys.exit("can not open outputfile")

        out.write(header+"\n")
        lines.remove(lines[0])

        for line in lines:
            line = line.rstrip("\n")
            splitted_line = line.split("\t")
            if add_column_for_novelty:
                splitted_line.append("")
            new_ncbiid = splitted_line[id_index]
            if new_ncbiid is "":
                continue
            new_ncbiid = int(new_ncbiid)
            novelty = self.get_novelty(new_ncbiid)

            new_line = ""
            if novelty is not None:
                print "\t" + str(new_ncbiid) + " is not included at rank " + str(novelty)
                splitted_line[novelty_index] = "new_" + str(novelty)

                for value in splitted_line:
                    new_line = new_line + value + "\t"
                new_line = new_line.rstrip("\t")
            else:
                new_line = line
            out.write(new_line)
            out.write("\n")
        out.close()

    def get_novelty(self,ncbiid):
        """
            Compute novelty_category for a new taxon ID according to the following algorithm
            @param ncbiid: NCBI ID for which the novelty_category should be computed
            @return: novelty category or None
            @rtype: str
        """

        if not self.tax.exists(ncbiid):
            print "\t[warning] " + str(ncbiid) + " not included in the taxonomy."
            return None

        taxons_rank = self.tax.getRank(ncbiid)
        try:
            taxons_rank_index = self.ranks.index(taxons_rank)
        except: #it is a strain
            if ncbiid in self.included_parents_at_rank['no rank']:
                return None     # may not even be a new strain
            else:
                return 'strain'

        for index in range(taxons_rank_index,len(self.ranks)-1):
            taxons_parent_at_rank = self.tax.parentAtRank(ncbiid,self.ranks[index])
            if taxons_parent_at_rank is not None:
                if taxons_parent_at_rank in self.included_parents_at_rank[self.ranks[index]]:
                    return self.ranks[index-1]

        return 'superkingdom'


    def get_taxonomic_ids_from_directory(self,dir):
        """
            search a directory for all files with taxonomic IDS and save them as a set
            @param dir: directory containing sequences named [ID].[nr].fna
            @return: set of the IDs
            @rtype: set
        """
        files = glob.glob(dir+"*.*")
        taxIDs = set()

        for f in files:
            tid = f.split("/")[-1].split(".")[0]
            if '_' in tid:
                tid = tid.spli('_')[0]
            taxIDs.add(int(tid))

        return taxIDs

    def close(self):
        self.tax.close()


#***************** Main ************************
if __name__ =='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-ref',action='store',help='dir with reference sequence files (named taxID.1.fna)',required=True)
    parser.add_argument('-db',action='store',help='complete path to ncbi-database file',required=True)
    parser.add_argument('-meta',action='store',required=True,help='metadatafile file')
    parser.add_argument('-o',action='store',required=True,help='output file')
    args = parser.parse_args()

    Nov = Novelty(args.db)
    Nov.read_reference(args.ref)
    Nov.compute_novelty_for_metafile(args.meta,args.o)
    Nov.close()
    print "Done."