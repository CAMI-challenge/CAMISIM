import gffutils
import random
import csv
import os
import argparse
import math

class GeneAbundace() :

    total_gene_size = 0

    def parse_annotationfile(self, db, number_of_samples, feature_type):
        """
            Get initial list with zero initialized

            @attention: Each list in the list contains the distribution value for each sample

            @param size_of_population: Amount of genomes or individuals
            @type size_of_population: int | long
            @param number_of_samples: Number of samples
            @type number_of_samples: int | long

            @return: A list of lists.
            @rtype: list[list[float]]
        """

        gene_id_to_abundances = dict()

        for gene in db.features_of_type(feature_type):
            # print(gene)
            gene_id = gene.id
            gene_id_to_abundances[gene_id] = [0.0] * number_of_samples

        return gene_id_to_abundances

    def add_initial_log_distribution(self, gene_id_to_abundances, mu, sigma, sample):
        """
            Adding a initial distribution

            @attention: Values for first sample

            @param gene_id_to_abundances: Main list for all distributions
            @type : list[list[float]]
            @param mu: Mean
            @type mu: float
            @param sigma: standard deviation
            @type sigma: float

            @return: Nothing
            @rtype: None
        """
        assert isinstance(gene_id_to_abundances, dict)
        assert isinstance(mu, (float, int))
        assert isinstance(sigma, (float, int))

        # ToDo
        # In Python 3.7.0 the insertion-order preservation nature of dict objects has been declared to be an official part of the Python language spec.
        # Therefore, you can depend on it.
        for key in gene_id_to_abundances.keys():
            gene_id_to_abundances[key][sample] = random.lognormvariate(mu, sigma)

        return gene_id_to_abundances

    def add_genewise_log_distribution(self, gene_id_to_abundances, sigma):

        assert isinstance(gene_id_to_abundances, dict)
        assert isinstance(sigma, (float, int))
        #ToDo wie genau distributed
    
        for gene in gene_id_to_abundances.keys():
            for sample in range(len(gene_id_to_abundances[gene])-1):
                mu_current_gene = gene_id_to_abundances[gene][sample]

                ln_mu_current_gene = math.log(mu_current_gene)
                gene_id_to_abundances[gene][sample+1] = random.lognormvariate(ln_mu_current_gene, sigma)

        return gene_id_to_abundances    

    def add_timeseries_gauss(self, gene_id_to_abundances, mu, sigma):
        """
            Adding gaussian noise sequentially to the previous sample

            @attention:

            @param gene_id_to_abundances: Main list for all distributions
            @type : list[list[float]]
            @param mu: Mean
            @type mu: float
            @param sigma: standard deviation
            @type sigma: float

            @return: Nothing
            @rtype: None
        """
        assert isinstance(gene_id_to_abundances, dict)
        assert isinstance(mu, (float, int))
        assert isinstance(sigma, (float, int))
        count = 0
        for key in gene_id_to_abundances.keys():
            for index_i in range(len(gene_id_to_abundances[key])-1):
                if gene_id_to_abundances[key][index_i] > 0:
                    gene_id_to_abundances[key][index_i+1] = self.lt_zero(gene_id_to_abundances[key][index_i] + random.gauss(mu, sigma))
                else:
                    # extinction
                    gene_id_to_abundances[key][index_i+1] = 0.0
                    count += 1
                    print("In extinction")
        print("Count of extinction: "+str(count))

        return gene_id_to_abundances


    def get_gene_size(self, db, gene_id):

        # Fetch the gene by ID
        gene = db[gene_id]

        # Calculate the size of the gene
        gene_size = gene.end - gene.start + 1

        return gene_size    

    def normalize_by_gene_size(self, db, gene_id_to_abundances):

        self.total_gene_size = 0

        for gene_id in list(gene_id_to_abundances.keys()):
            size = self.get_gene_size(db, gene_id)

            self.total_gene_size = self.total_gene_size + size

            for i in range(len(gene_id_to_abundances[gene_id])):

                gene_id_to_abundances[gene_id][i] *= size
                
                # Test für nanosin
                # gene_id_to_abundances[gene_id][i] /= (size/1000)

        return gene_id_to_abundances        

    def random_distribution_to_relative_abundance(self, gene_id_to_abundances, number_of_samples, precision=10):
        """
            Replace random distributions with relative abundances

            @attention: limited to first 20 ???

            @param list_population: Main list for all distributions
            @type list_population: list[list[float]]
            @param precision: Precision, numbers after decimal point
            @type precision: int
        """

        for i in range(number_of_samples):

            total_abundance = 0.0

            for gene_id in list(gene_id_to_abundances.keys()):
                total_abundance += gene_id_to_abundances[gene_id][i]

            for gene_id in list(gene_id_to_abundances.keys()):
                gene_id_to_abundances[gene_id][i] = round(gene_id_to_abundances[gene_id][i] / float(total_abundance), precision)

        return gene_id_to_abundances               

    def normalize_abundances(self, db, gene_id_to_abundances, number_of_samples):

        gene_id_to_abundances = self.normalize_by_gene_size(db, gene_id_to_abundances)

        gene_id_to_abundances = self.random_distribution_to_relative_abundance(gene_id_to_abundances, number_of_samples)

        return gene_id_to_abundances

    def write_dict_to_tsv(self, data, file_basename, number_of_samples, precision=10):

        for i in range(number_of_samples):

            filename = file_basename + "_sample_" + str(i) + ".tsv"

            # Open the file in write mode
            with open(filename, 'w', newline='') as file:
                # Iterate over the dictionary items
                for key, value in data.items():
                    # Format the value as a float with variable precision
                    formatted_value = f"{float(value[i]):.{precision}f}"
                    # Write the key and formatted value separated by a tab, then a newline
                    file.write(f"{key}\t{formatted_value}\n")

    def lt_zero(self, value):
        """
            Prevent values <= 0

            @attention:

            @param value:
            @type value: float | int | long

            @return: value if > 0, else 0.001
            @rtype: float | int | long
        """
        if value <= 0:
            # > 0 to prevent extinction
            return 0.001
        else:
            return value

    def create_db(self, annotation_file):

        return gffutils.create_db(annotation_file, dbfn=(annotation_file+".db"), merge_strategy="create_unique", keep_order=True)

    def _get_gene_abundance(self):

        parser = argparse.ArgumentParser()
        parser.add_argument(
		    "-annotation_file",
		    help="the annotation file in gff3 format",
		    action='store',
            required=True)
        parser.add_argument(
		    "-mode",
		    help="the mode of the distribution",
		    action='store',
		    required=True)
        parser.add_argument(
            "-number_of_samples",
            help="the number of samples",
            type=int,
            required=True)
        parser.add_argument(
            "-seed",
            help="the seed to ensure reproducibility",
            type=int,
            required=True)
        parser.add_argument(
            "-log_mu",
            help="The mean to use for the lognormal distribution",
            type=float,
            required=True)
        parser.add_argument(
            "-log_sigma",
            help="The standard deviation to use for the lognormal distribution",
            type=float,
            required=True)
        parser.add_argument(
            "-gauss_mu",
            help="The mean to use for the gaussian noise",
            type=float,
            required=False)
        parser.add_argument(
            "-gauss_sigma",
            help="The standard deviation to use for the gaussian noise",
            type=float,
            required=False)
        parser.add_argument(
            "-gene_sigma",
            help="The standard deviation to use for the gaussian noise",
            type=float,
            required=False)
        parser.add_argument(
            '-feature_type',
            action='store',
            help='A feature type',
            required=True)       
        options = parser.parse_args()

        annotation_file = options.annotation_file
        mode = options.mode
        number_of_samples = options.number_of_samples
        seed = options.seed
        log_mu = options.log_mu
        log_sigma = options.log_sigma
        feature_type = options.feature_type

        if seed is not None:
                random.seed(seed)

        db = self.create_db(annotation_file)

        gene_id_to_abundances = self.parse_annotationfile(db, number_of_samples, feature_type)

        gene_id_to_abundances = self.add_initial_log_distribution(gene_id_to_abundances, log_mu, log_sigma, sample=0)

        if(mode == "differential"):

            for sample in range(1, number_of_samples):

                gene_id_to_abundances = self.add_initial_log_distribution(gene_id_to_abundances, log_mu, log_sigma, sample=sample)

        elif(mode == "timeseries"):

            gauss_mu = options.gauss_mu
            gauss_sigma = options.gauss_sigma

            gene_id_to_abundances = self.add_timeseries_gauss(gene_id_to_abundances, gauss_mu, gauss_sigma)

        elif(mode == "replicate"):

            gene_sigma = options.gene_sigma

            gene_id_to_abundances = self.add_genewise_log_distribution(gene_id_to_abundances, gene_sigma)
            # gene_id_to_abundances = self.add_genewise_log_distribution(self.normalize_abundances(db, gene_id_to_abundances, 1), gene_sigma)
        

        gene_id_to_abundances = self.normalize_abundances(db, gene_id_to_abundances, number_of_samples)

        base_file_name = "distribution_" + os.path.basename(annotation_file)

        self.write_dict_to_tsv(gene_id_to_abundances, base_file_name, number_of_samples)

        print(str(self.total_gene_size))


# main method and entry point of this script
if __name__ == "__main__":

    GeneAbundace()._get_gene_abundance()

