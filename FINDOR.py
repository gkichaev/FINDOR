
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.stats import chi2
import argparse


class FINDOR:

    def __init__(self, args):
        self.args = args
        self.taus_estimates = pd.Series()
        self.gwas_data = pd.DataFrame()
        self.predicted_tag_var = pd.DataFrame(columns=['SNP', 'PTV'])

    def tagged_variance_single(self, curr_fname):
        try:
            ld_scores = pd.read_table(curr_fname, header=0, delim_whitespace=True)
            pred_tagged_var = pd.DataFrame(columns=['SNP', 'PTV'])
            pred_tagged_var['SNP'] = ld_scores['SNP']
            just_ld_scores = ld_scores.drop(['CHR', 'SNP', 'BP'], axis=1)
            pred_tagged_var['PTV'] = np.dot(just_ld_scores,self.taus_estimates)
            del ld_scores
            self.predicted_tag_var = pd.concat([self.predicted_tag_var, pred_tagged_var])
        except:
            print("WARNING: LD-score file", curr_fname, "does not exists.",
                  "\nSNPs corresponding to this file will be excluded from analysis")

    def read_ldscores_chr(self, fname_prefix,  gzipped=True):
        if(gzipped):
            suffix = ".l2.ldscore.gz"
        else:
            suffix= ".l2.ldscore"
        names = {'prefix':fname_prefix, 'suffix':suffix}
        print('Reading in chromosome specific LD-scores from files: {prefix}[1-22]{suffix}'.format(**names))
        for i in range(1,23):
            curr_fname = fname_prefix+str(i)+suffix
            self.tagged_variance_single(curr_fname)

    def create_bins(self, num_bins):
        self.intersect_SNPs()
        self.gwas_data['PTV'] =  self.gwas_data['PTV']*self.gwas_data['N']
        quantiles = np.arange(0, 1 + 1 / num_bins, 1 / num_bins)
        predicted_tagged_variance_quantiles = self.gwas_data['PTV'].quantile(quantiles)
        #expand top quantiles to ensure everything is within range
        predicted_tagged_variance_quantiles[0] = predicted_tagged_variance_quantiles[0]-1
        predicted_tagged_variance_quantiles[1] = predicted_tagged_variance_quantiles[1]+1

        bins = pd.cut(self.gwas_data['PTV'], predicted_tagged_variance_quantiles, labels=np.arange(num_bins)) #create the lables
        self.gwas_data['bin_number'] = bins
        self.predicted_tag_var = None #no longer need this object

    def storey_pi_estimator(self, bin_index):
        """
        Estimate pi0/pi1 using Storey and Tibshirani (PNAS 2003) estimator.
        Argss
        =====
        bin_index: array of indices for a particular bin
        """
        pvalue = self.gwas_data.loc[bin_index,'P'] # extract pvalues from specific bin based index
        assert(pvalue.min() >= 0 and pvalue.max() <= 1), "Error: p-values should be between 0 and 1"
        total_tests = float(len(pvalue))
        pi0 = []
        lam = np.arange(0.05, 0.95, 0.05)
        counts = np.array([(pvalue > i).sum() for i in np.arange(0.05, 0.95, 0.05)])
        for l in range(len(lam)):
            pi0.append(counts[l] / (total_tests * (1 - lam[l])))

        # fit  cubic spline
        cubic_spline = interpolate.CubicSpline(lam, pi0)
        pi0_est = cubic_spline(lam[-1])
        if(pi0_est >1): #take care of out of bounds estimate
            pi0_est = 1
        return pi0_est

    def reweight_pvalue(self,num_bins):
        self.gwas_data['pi0'] = None
        print("Estimating pi0 within each bin")
        for i in range(num_bins):
            bin_index = self.gwas_data['bin_number']== i # determine index of snps in bin number i
            pi0 = self.storey_pi_estimator(bin_index)
            self.gwas_data.loc[bin_index, 'pi0'] = pi0
        if any(self.gwas_data['pi0'] == 1): # if a bin is estimated to be all null, give the smallest non-null weight
            one_index = self.gwas_data['pi0'] == 1
            largest_pi0 = self.gwas_data.loc[~one_index]['pi0'].max()
            self.gwas_data.loc[one_index,'pi0'] = largest_pi0
        print("Re-weighting SNPs")
        weights = (1-self.gwas_data['pi0'])/(self.gwas_data['pi0'])
        mean_weight = weights.mean()
        weights = weights/mean_weight #normalize weights to have mean 1
        self.gwas_data['weights'] = weights
        self.gwas_data['P_weighted'] = self.gwas_data['P']/weights #reweight SNPs

    def read_sumstats(self, fname):
        print("Reading in GWAS data from file", fname)
        self.gwas_data = pd.read_table(fname, header = 0, delim_whitespace=True)
        missing = [i not in self.gwas_data for i in ["SNP", "Z", "N"]]
        if any(missing):
            raise ValueError('Need SNP, Z, and N fields in the GWAS data file')
        self.gwas_data= self.gwas_data.dropna(axis=0, how='any')

    def zscores_to_pvalues(self):
        if('P' not in self.gwas_data.columns):
            chi_squares = self.gwas_data['Z'].apply(np.square)
            pvalues = chi2.sf(chi_squares, df=1)
            self.gwas_data['P'] = pvalues
        else:
            pass

    def read_taus(self, fname):
        print("Reading LD-score regression results from file", fname)
        full_taus = pd.read_table(fname, header = 0, delim_whitespace=True)
        self.taus_estimates = full_taus['Coefficient']

    def write_output(self, fname):
        print('Writing output to:', fname)
        index = self.gwas_data['P_weighted'] > 1
        self.gwas_data.loc[index, 'P_weighted'] = 1
        if(self.args.bin_info):
            self.gwas_data.to_csv(fname, sep = ' ', index=False)
        else:
            self.gwas_data.drop(['PTV', 'bin_number' ,'pi0' ,'weights'], axis=1).to_csv(fname, sep=' ', index=False)

    def intersect_SNPs(self):
        """
        Intersect SNPS with LD scores and (optionally) with a SNP list provided
        :param:
        :return:
        """
        print("Intersecting GWAS SNPs with LD scores")
        self.gwas_data = pd.merge(self.gwas_data, self.predicted_tag_var, on='SNP')
        print("After intersection with LD scores", str(self.gwas_data.shape[0]), "remain")
        fname = self.args.snp_list
        if(fname is not None):
            print("Reading in SNP list from file", fname)
            snp_list = pd.read_table(fname, header = 0, delim_whitespace=True)
            if('SNP' in snp_list.columns):
                self.gwas_data = pd.merge(self.gwas_data, snp_list, on='SNP')
                num_remain = self.gwas_data.shape
                print("After merging with SNP list", str(num_remain[0]), "SNPs remain in data")
            else:
                ValueError("SNP list does not header with 'SNP'")

    def run_algorithm(self):

        if (self.args.ldscore_path != None and self.args.ldscore_path_chr == None):
            path_to_ldscores = self.args.ldscore_path
        elif (self.args.ldscore_path == None and self.args.ldscore_path_chr != None):
            path_to_ldscores = self.args.ldscore_path_chr
        else:
            raise ValueError("Error: Need to specify path to LD scores")

        self.read_taus(self.args.tau_path)
        if(path_to_ldscores[-1] == "."):
            self.read_ldscores_chr(path_to_ldscores)
        else:
            self.tagged_variance_single(path_to_ldscores)
        num_ldscores = self.predicted_tag_var.shape
        print("Read in LD-scores for" , str(num_ldscores[0]), "SNPs")
        self.read_sumstats(self.args.gwas_path)
        num_gwas_snps = self.gwas_data.shape
        print("Read in GWAS data for", str(num_gwas_snps[0]), "SNPs")
        self.zscores_to_pvalues()
        self.create_bins(self.args.num_bins)
        self.reweight_pvalue(num_bins = self.args.num_bins)
        self.write_output(self.args.out_file)

# get command line
def parse_command_line():
    parser = argparse.ArgumentParser(description='Functionally Informed Novel  Discovery Of Risk Loci (FINDOR)\n'
                                                 'Gleb Kichaev 2017\n'
                                                 'Please send bug-reports to gkichaev@ucla.edu')

    parser.add_argument('--ref-ld', default=None, type=str, dest='ldscore_path',
                        help='Use --ref-ld to specify the path to the LD scores file to use to predict the tagged variance.'
                             ' Note: the dimensions of this file must align with the ".results" file. '
                             )
    parser.add_argument('--ref-ld-chr', default=None, type=str, dest="ldscore_path_chr",
                        help='Use --ref-ld-chr to specify the prefix to the LD score files that are separated by chromosome. For example,'
                             'if LD scores for chromsomes 1-22 are stored in LD/baselineLD.{1..22}.l2.ldscore.gz, then use'
                             '"--ref-ld-chr LD/baselineLD."'
                             'Note that the suffix, l2.ldscore.gz is automatically appended\n.')
    parser.add_argument('--gwas-data',default=None, type=str, dest='gwas_path',
                        help = '(required) Use --gwas-data to specify the path to GWAS data file (must have SNP, Z, N headers).\n',
                        required=True)
    parser.add_argument('--regression-results', default=None, type=str, dest='tau_path',
                        help= '(required) use --regression-results to specify the path to the .results file'
                              'from applying LD-score regression on your data.'
                        'These values contain the coefficients for effect sizes of the annotations.\n', required=True)
    parser.add_argument('--num-bins', default=100, type=int, dest='num_bins',
                        help = 'Use the --num-bins flag to specify the number of bins to use when stratifying SNPs (default = 100)\n')
    parser.add_argument('--out', dest='out_file', type=str,
                        help='(required) use --out to specify the output file name for re-weighted\n', required=False)
    parser.add_argument('--output_bin_info', dest='bin_info', default=False, type=bool,
                        help='Should bin information (PTV, bin_number, pi0, weights) be included in the output? (Default=False)')
    parser.add_argument('--snp-list', dest='snp_list', default=None, type=str,
                        help='Path to a list of SNPs to analyze ("SNP" header expected)? (Default=None)')
    args = parser.parse_args()
    return args


def main():
    args = parse_command_line()
    findor = FINDOR(args)
    findor.run_algorithm()

if __name__ == "__main__":
    main()