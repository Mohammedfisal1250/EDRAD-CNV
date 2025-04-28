from cbs import segment
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from PIL import ImageTk, Image
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import pysam
import pandas as pd
import datetime
import math
import sys
from sklearn.cluster import KMeans
import scipy
import gc
import os.path
from sklearn.preprocessing import StandardScaler
from scipy.stats import entropy
from sklearn.neighbors import KernelDensity
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import NearestNeighbors
import rpy2.robjects as robjects
from scipy import stats
import warnings

# [All your existing functions remain exactly the same...]
# Define functions here (the same functions from your code)

def run_analysis(bam, ref_path, binSize, outfile, p_value_file, col, k_value, groudtruth, result_file, score_result_file, segmentation_method, bandwidth):
    try:
    # Your main code here
        def calculate_density_gmm(X, n_components=1):
            gmm = GaussianMixture(n_components=n_components)
            gmm.fit(X)
            log_density = gmm.score_samples(X)
            return np.exp(log_density)

        def calculate_density_knn(X, k=5):
            nbrs = NearestNeighbors(n_neighbors=k).fit(X)
            distances, _ = nbrs.kneighbors(X)
            densities = 1.0 / (np.mean(distances, axis=1) + 1e-10)  # Inverse of average distance to k-nearest neighbors
            return densities

        def calculate_density_histogram(X, bins=10):
            hist, bin_edges = np.histogram(X, bins=bins, density=True)
            bin_indices = np.digitize(X.flatten(), bin_edges) - 1
            densities = hist[bin_indices]
            return densities

        def calculate_density(X, bandwidth):
            kde = KernelDensity(bandwidth=bandwidth)
            kde.fit(X)
            log_density = kde.score_samples(X)
            return np.exp(log_density)

        def calculate_local_entropy(densities, k):
            local_entropy = np.zeros(densities.shape)
            for i, density in enumerate(densities):
                distances = np.abs(densities - density)
                k_nearest_densities = np.partition(distances, k)[:k]
                local_entropy[i] = entropy(k_nearest_densities + 1e-10)  # Add small value to avoid log(0)
            return local_entropy

        def edrad(X, k,  method='kde', bandwidth=bandwidth, n_components=1):
            if method == 'kde':
                densities = calculate_density(X, bandwidth)
            elif method == 'gmm':
                densities = calculate_density_gmm(X, n_components)
            elif method == 'knn':
                densities = calculate_density_knn(X, k)
            elif method == 'histogram':
                densities = calculate_density_histogram(X)
            else:
                raise ValueError("Unknown method")
            
            local_entropies = calculate_local_entropy(densities, k)
            edr_scores = local_entropies / (densities + 1e-10)  # Add small value to avoid division by zero
            return edr_scores


        def normalize(scores):
            return (scores - np.min(scores)) / (np.max(scores) - np.min(scores))

        def calc_scores(data,k, bandwidth, method='kde'):
            # Calculate EDR scores using edrad
            edr_scores = edrad(data, k=k, method=method, bandwidth=bandwidth, n_components=1)  # Adjust parameters based on your requirements
            
            # Normalize scores
            normalized_scores = normalize(edr_scores)
            
            return normalized_scores
        
        def read_bam(file):
            try:
                # reading bam file
                print(f"Opening BAM file: {file}")
                sam_file = pysam.AlignmentFile(file, "rb", ignore_truncation=True)
                chr_list = sam_file.references
                print(f"Chromosomes in BAM file: {chr_list}")
                return chr_list
            except (IOError, ValueError, OSError) as e:
                print(f"Error reading BAM file: {e}")
                raise

        def read_ref(file, chr_num, ref):
            # read reference file
            if os.path.exists(file):
                print("Read reference file: " + str(file))
                with open(file, 'r') as f:
                    line = f.readline()
                    for line in f:

                        lines = line.strip()
                        ref[chr_num] += lines
            else:
                print("Warning: can not open " + str(file) + '\n')
            return ref

        def bins(ref, bin_size, chr_len, file):
            chr_tag = np.full(23, 0)
            chr_list = np.arange(23)
            max_num = int(chr_len.max() / bin_size) + 1
            init_rd = np.full((23, max_num), 0.0)
            try:
                # read bam file and get bin rd
                print("Read bam file: " + str(file))
                sam_file = pysam.AlignmentFile(file, "rb", ignore_truncation=True)
                for line in sam_file:
                    idx = int(line.pos / bin_size)
                    if idx > int(chr_len.max() // bin_size):
                        continue
                    if line.reference_name:
                        chr = line.reference_name.strip('chr')
                        if chr.isdigit():
                            init_rd[int(chr)][idx] += 1
                            chr_tag[int(chr)] = 1
            except (IOError, ValueError, OSError) as e:
                print(f"Error reading BAM file: {e}")
                raise
            except pysam.utils.PysamError as e:
                print(f"Pysam error: {e}")
                raise

            chr_list = chr_list[chr_tag > 0]
            chr_num = len(chr_list)
            rd_list = [[] for _ in range(chr_num)]
            pos_list = [[] for _ in range(chr_num)]
            init_gc = np.full((chr_num, max_num), 0)
            pos = np.full((chr_num, max_num), 0)
            # initialize bin_data and bin_head
            count = 0
            for i in range(len(chr_list)):
                chr = chr_list[i]
                bin_num = int(chr_len[chr] / bin_size) + 1
                for j in range(bin_num):
                    pos[i][j] = j
                    cur_ref = ref[chr][j * bin_size:(j + 1) * bin_size]
                    N_count = cur_ref.count('N') + cur_ref.count('n')
                    if N_count == 0:
                        gc_count = cur_ref.count('C') + cur_ref.count('c') + cur_ref.count('G') + cur_ref.count('g')
                    else:
                        gc_count = 0
                        init_rd[chr][j] = -1000000
                        count = count + 1
                    init_gc[i][j] = int(round(gc_count / bin_size, 3) * 1000)
                # delete
                cur_rd = init_rd[chr][:bin_num]
                cur_gc = init_gc[i][:bin_num]
                cur_pos = pos[i][:bin_num]
                cur_rd = cur_rd / 1000
                index = cur_rd >= 0
                rd = cur_rd[index]
                GC = cur_gc[index]
                cur_pos = cur_pos[index]
                rd[rd == 0] = mode_rd(rd)
                rd = gc_correct(rd, GC)
                pos_list[i].append(cur_pos)
                rd_list[i].append(rd)
            del init_rd, init_gc, pos
            gc.collect()
            return rd_list, pos_list, chr_list


        def mode_rd(rd):
            new_rd = np.full(len(rd), 0)
            for i in range(len(rd)):
                new_rd[i] = int(round(rd[i], 3) * 1000)
            count = np.bincount(new_rd)
            count_list = np.full(len(count) - 49, 0)
            for i in range(len(count_list)):
                count_list[i] = np.mean(count[i:i + 50])
            mode_min = np.argmax(count_list)
            mode_max = mode_min + 50
            mode = (mode_max + mode_min) / 2
            mode = mode / 1000
            return mode


        def gc_correct(rd, gc):
            # correcting gc bias
            bin_count = np.bincount(gc)
            global_rd_ave = np.mean(rd)
            for i in range(len(rd)):
                if bin_count[gc[i]] < 2:
                    continue
                mean = np.mean(rd[gc == gc[i]])
                rd[i] = global_rd_ave * rd[i] / mean
            return rd

        def scaling_rd(rd, mode):
            posit_rd = rd[rd > mode]
            neg_rd = rd[rd < mode]
            if len(posit_rd) < 50:
                mean_max_rd = np.mean(posit_rd)
            else:
                sort = np.argsort(posit_rd)
                max_rd = posit_rd[sort[-50:]]
                mean_max_rd = np.mean(max_rd)
            if len(neg_rd) < 50:
                mean_min_rd = np.mean(neg_rd)
            else:
                sort = np.argsort(neg_rd)
                min_rd = neg_rd[sort[:50]]
                mean_min_rd = np.mean(min_rd)
            scaling = mean_max_rd / (mode + mode - mean_min_rd)
            for i in range(len(rd)):
                if rd[i] < mode:
                    rd[i] /= scaling
            return rd


        def plot(pos, data):
            pos1 = np.arange(1, len(data)+1)
            plt.scatter(pos1, data, s=5, c="black")
            plt.xlabel("pos")
            plt.ylabel("scalRD")
            plt.show()

        def plotgc(pos, data):
            plt.scatter(pos, data, s=3, c="blue")
            #plt.scatter(pos1, data1, s=3, c="red")
            plt.xlabel("ds")
            plt.ylabel("gc")
            plt.show()

        def plotGCscatter(gc,rd):
            plt.scatter(gc, rd, s=1,linewidths=0.1)
            plt.xlabel("GC_content")
            plt.ylabel("rd")
            plt.show()

        def seg_rd(rd, bin_head, seg_start, seg_end, seg_count):
            seg_rd = np.full(len(seg_count), 0.0)
            for i in range(len(seg_rd)):
                seg_rd[i] = np.mean(rd[seg_start[i]:seg_end[i]])
                seg_start[i] = bin_head[seg_start[i]] * binSize + 1
                if seg_end[i] == len(bin_head):
                    seg_end[i] = len(bin_head) - 1
                seg_end[i] = bin_head[seg_end[i]] * binSize + binSize
            return seg_rd, seg_start, seg_end
        
        def segmentation_cbs_py(rd, pos, binSize):
            #messagebox.showinfo("Selection", "Python script selected: segmentation_cbs_py")
            def _get_rd_values(rd, pos, seg_index, binSize):
                 
                seg_rd = []
                seg_start = np.full(len(seg_index), 0)
                seg_end = np.full(len(seg_index), 0)
                for i in range(len(seg_index)):
                    segment = rd[seg_index[i][0]: seg_index[i][1]]
                    seg_rd.append([np.mean(segment)])
                    seg_start[i] = pos[seg_index[i][0]] * binSize + 1
                    if seg_end[i] == len(pos):
                       seg_end[i] = len(pos) - 1
                    seg_end[i] = pos[seg_index[i][1] - 1] * binSize + binSize

                return seg_rd, seg_start, seg_end

            seg_index = segment(scale_rd)
            return _get_rd_values(rd, pos, seg_index, binSize)
        
        def segmentation_cbs_r(seg_path, rd, pos, binSize, num_bin, col=col):
            #messagebox.showinfo("Selection", "R script selected: segmentation_cbs_r")
            def _get_rd_values(rd, pos, seg_start, seg_end, binSize):
                per_seg_rd = []
                for i in range(len(seg_end)):
                    seg = rd[seg_start[i]: seg_end[i]]
                    per_seg_rd.append([np.mean(seg)])
                    seg_start[i] = pos[seg_start[i]] * binSize + 1
                    if seg_end[i] == len(pos):
                        seg_end[i] = len(pos) - 1
                    seg_end[i] = pos[seg_end[i]] * binSize + binSize

                return per_seg_rd, seg_start, seg_end

            v = robjects.FloatVector(scale_rd)
            m = robjects.r['matrix'](v, ncol=col)
            robjects.r.source("CBS_data.R")
            robjects.r.CBS_data(m, seg_path)

            num_col = int(num_bin / col) + 1
            seg_start, seg_end, seg_count, seg_len = read_seg_file(num_col, num_bin)

            seg_start = seg_start[:-1]
            seg_end = seg_end[:-1]

            return _get_rd_values(rd, pos, seg_start, seg_end, binSize)

        def write_data_file(chr, seg_start, seg_end, seg_count, outlier_scores,label,index):
            output = open(p_value_file, "w")
            output.write("Chr Num " + '\t' + " Start Position " + '\t' + " End Position " + '\t' + "  all_rd " + '\t\t' + " outlier_scores " + '\t\t' +"label"+'\t\t' +"index"+'\n')

            for i in range(len(outlier_scores)):
                output.write(str(chr[i]) + '\t ' + str(seg_start[i]) + ' \t ' + str(seg_end[i]) + ' \t ' + str(seg_count[i]) + ' \t ' + str(outlier_scores[i]) + ' \t ' +str(label[i])+' \t ' +str(index[i])+ '\n')

        def write_dataset_file(seg_start, seg_count,label):
            """
            write dataset file
            outlier_scores,label
            """
            output_dataset = open(dataset_file, "w")
            output_dataset.write("seg_start " + '\t' + " all_rd " + '\t' +"label"+'\n')

            for i in range(len(outlier_scores)):
                output_dataset.write( str(seg_start[i]) + '\t ' + str(seg_count[i]) + ' \t ' +str(label[i])+ '\n')


        def write_cnv_file(chr, cnv_start, cnv_end, cnv_type, cn, filename):
            """
            write cnv result file
            pos start, pos end, type, copy number
            """
            output = open(filename, "w")
            for i in range(len(cnv_type)):
                if cnv_type[i] == 2:
                    output.write("Chr" + str(chr[i]) + '\t' + str(cnv_start[i]) + '\t' + str(
                        cnv_end[i]) + '\t' + str("duplication") + '\t' + str(cn[i]) + '\n')
                else:
                    output.write("Chr" + str(chr[i]) + '\t' + str(cnv_start[i]) + '\t' + str(
                        cnv_end[i]) + '\t' + str("deletion") + '\t' + str(cn[i]) + '\n')


        def read_seg_file(num_col, num_bin):
            """
            read segment file (Generated by DNAcopy.segment)
            seg file: col, chr, start, end, num_mark, seg_mean
            """
            seg_start = []
            seg_end = []
            seg_count = []
            seg_len = []
            with open("seg", 'r') as f:
                for line in f:
                    line_str_list = line.strip().split('\t')
                    start = (int(line_str_list[0]) - 1) * num_col + int(line_str_list[2]) - 1
                    end = (int(line_str_list[0]) - 1) * num_col + int(line_str_list[3]) - 1
                    if start < num_bin:
                        if end > num_bin:
                            end = num_bin - 1
                        seg_start.append(start)
                        seg_end.append(end)
                        seg_count.append(float(line_str_list[5]))
                        seg_len.append(int(line_str_list[4]))
            seg_start = np.array(seg_start)
            seg_end = np.array(seg_end)
            return seg_start, seg_end, seg_count, seg_len


        def calculating_copy_number(mode, cnv_rd, cnv_type):
            cn = np.full(len(cnv_type), 0)
            index = cnv_type == 1
            lossRD = cnv_rd[index]
            if len(lossRD) > 2:
                data = np.c_[lossRD, lossRD]
                del_type = KMeans(n_clusters=2, random_state=9).fit_predict(data)
                cnv_type[index] = del_type
                if np.mean(lossRD[del_type == 0]) < np.mean(lossRD[del_type == 1]):
                    homo_rd = np.mean(lossRD[del_type == 0])
                    hemi_rd = np.mean(lossRD[del_type == 1])
                    for i in range(len(cn)):
                        if cnv_type[i] == 0:
                            cn[i] = 0
                        elif cnv_type[i] == 1:
                            cn[i] = 1
                else:
                    hemi_rd = np.mean(lossRD[del_type == 0])
                    homo_rd = np.mean(lossRD[del_type == 1])
                    for i in range(len(cn)):
                        if cnv_type[i] == 1:
                            cn[i] = 0
                        elif cnv_type[i] == 0:
                            cn[i] = 1
                purity = 2 * (homo_rd - hemi_rd) / (homo_rd - 2 * hemi_rd)
                for i in range(len(cnv_type)):
                    if cnv_type[i] == 2:
                        # Extract scalar value from cnv_rd[i] if it is an array
                        cnv_rd_i = cnv_rd[i] if np.isscalar(cnv_rd[i]) else cnv_rd[i].item()
                        cn[i] = int(2 * cnv_rd_i / (mode * purity) - 2 * (1 - purity) / purity)
            return cn


        def boxplot(outlier_scores):
            four = pd.Series(outlier_scores).describe()
            Q1 = four['25%']
            Q3 = four['75%']
            IQR = Q3 - Q1
            upper = Q3 + 0.75 * IQR
            lower = Q1 - 0.75 * IQR
            return upper

        def combining_cnv(seg_chr, seg_start, seg_end, seg_count,outlier_scores, upper, mode):
            index = outlier_scores > upper
           # index = outlier_scores > 0.85
            print("index=",index)
            CNV_chr = seg_chr[index]
            CNV_start = seg_start[index]
            CNV_end = seg_end[index]
            CNV_RD = seg_count[index]
            #CNV_label = seg_len[index]
            type = np.full(len(CNV_RD), 1)
            for i in range(len(CNV_RD)):
                if CNV_RD[i] > mode:
                    type[i] = 2
            for i in range(len(CNV_RD) - 1):
                if CNV_end[i] + 1 == CNV_start[i + 1] and type[i] == type[i + 1]:
                    CNV_start[i + 1] = CNV_start[i]
                    type[i] = 0
            index = type != 0
            CNV_RD = CNV_RD[index]
            CNV_chr = CNV_chr[index]
            CNV_start = CNV_start[index]
            CNV_end = CNV_end[index]
            CNV_type = type[index]
            return CNV_chr, CNV_start, CNV_end, CNV_RD, CNV_type
        
        def label_score(outlier_scores, upper):
            label = np.full(len(outlier_scores), 0)
            print (label)
            for i in range(len(outlier_scores)):
                if outlier_scores[i] > upper:
                    label[i]=1
                else:
                    label[i]=0
            return label

        def calculate_scores(groudtruth, result_file, score_result_file):
            result_start = []
            result_end = []
            result_type = []

            with open(result_file, 'r') as f:
                for line in f:
                    linestr = line.strip()
                    linestrlist = linestr.split('\t')
                    result_start.append(int(linestrlist[1]))
                    result_end.append(int(linestrlist[2]))
                    result_type.append(linestrlist[3])

            ground_truth = pd.read_table(groudtruth)
            truth_type = ground_truth["variant type"].tolist()
            truth_start = ground_truth['start'].tolist()
            truth_end = ground_truth['stop'].tolist()

            count = 0
            for i in range(len(result_type)):
                for j in range(len(truth_type)):
                    if truth_start[j] <= result_start[i] <= truth_end[j] and truth_type[j] == result_type[i]:
                        if result_end[i] <= truth_end[j]:
                            count += (result_end[i] - result_start[i] + 1)
                        elif result_end[i] > truth_end[j]:
                            count += (truth_end[j] - result_start[i] + 1)
                    elif truth_start[j] >= result_start[i] and truth_type[j] == result_type[i]:
                        if truth_start[j] <= result_end[i] <= truth_end[j]:
                            count += (result_end[i] - truth_start[j] + 1)
                        elif result_end[i] >= truth_end[j]:
                            count += (truth_end[j] - truth_start[j] + 1)

            result_count = sum(result_end[i] - result_start[i] + 1 for i in range(len(result_start)))
            truth_count = sum(truth_end[i] - truth_start[i] + 1 for i in range(len(truth_start)))

            precision = count / result_count if result_count != 0 else 0
            sensitivity = count / truth_count if truth_count != 0 else 0
            F1_score = (2 * precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) != 0 else 0

            with open(score_result_file, 'w') as f4:
                f4.write('result_' + "sensitivity" + ' = '+ str(sensitivity) + '\n')
                f4.write('result_' + "precision" + ' = '+ str(precision) + '\n')
                f4.write('result_' + "F1_score" + ' = '+ str(F1_score) + '\n')

            return sensitivity, precision, F1_score

        print("BAM File:", bam)
        print("Reference Path:", ref_path)
        print("Bin Size:", binSize)
        print("Output File:", outfile)
        print("P-Value File:", p_value_file)
        print("Column:", col)
        #print("Dataset File:", dataset_file)
        #print("Data Filename:", data_filename)
        print("k_value:", k_value)
        print("bandwidth:", bandwidth)
        print("groudtruth_path", groudtruth)
        print("result_file", result_file)
        print("score_result_file", score_result_file)
        print("Segmentation Method:", segmentation_method)
        
        start = datetime.datetime.now()
        path = os.path.abspath('.')
        seg_path = path + str("/seg")

        ref = [[] for i in range(23)]
        refList = read_bam(bam)
        for i in range(len(refList)):
            chr = refList[i]
            if chr == '21':
                fa_seq = SeqIO.read(ref_path, "fasta")
                ref[21] = str(fa_seq.seq)

        chrLen = np.full(23, 0)
        for i in range(1, 23):
            chrLen[i] = len(ref[i])
        RDList, PosList, chrList = bins(ref, binSize, chrLen, bam)
        all_chr = []
        all_rd = []
        all_start = []
        all_end = []
        modeList = np.full(len(chrList), 0.0)
        for i in range(len(chrList)):
           # print("analyse " + str(chrList[i]))
            RD = np.array(RDList[i][0])
            pos = np.array(PosList[i][0])
            num_bin = len(RD)
            modeList[i] = mode_rd(RD)
            scale_rd = scaling_rd(RD, modeList[i])
            print("segment count...")

            # Call the appropriate segmentation function based on user selection
            if segmentation_method == "Python":
                seg_rd, seg_start, seg_end = segmentation_cbs_py(scale_rd, pos, binSize)
            else:
                seg_rd, seg_start, seg_end = segmentation_cbs_r(seg_path, scale_rd, pos, binSize, num_bin, col)

            all_rd.extend(seg_rd)
            all_start.extend(seg_start)
            all_end.extend(seg_end)
            all_chr.extend(chrList[i] for _ in range(len(seg_rd)))

        all_chr = np.array(all_chr)
        all_start = np.array(all_start)
        all_end = np.array(all_end)
        all_rd = np.array(all_rd)
        for i in range(len(all_rd)):
            if np.isnan(all_rd[i]).any():
                all_rd[i] = (all_rd[i - 1] + all_rd[i + 1]) / 2


        print("Calculating scores...")
        print("all_rd=", all_rd)
        data1 = all_rd.reshape(-1, 1)
        data = StandardScaler().fit_transform(data1)

        print("data size=", data.size)
        print("data=", data)

        #data_filename = "19239_edrad_GRCh38.dna_data.csv"
        #df = pd.DataFrame(data, columns=['RD_Scores'])  # Assuming 'data' is a 2D array with one column
        #df.to_csv(data_filename, index=False)

        #print(f"Data saved to {data_filename}")

        #outlier_scores = calc_scores(data)
        outlier_scores = calc_scores(data, k_value, bandwidth)
        #outlier_scores = calc_scores(data, method='gmm')
        #outlier_scores = calc_scores(data, method='knn')
        #outlier_scores = calc_scores(data, method='histogram')
        print("outlier_scores size =", outlier_scores.size)


        upper = boxplot(outlier_scores)
        #upper = Q3 + 0.75 * IQR
        label = label_score(outlier_scores, upper)
        print ("label=",label)

        index = outlier_scores > upper
        mode = np.mean(modeList)
        write_data_file(all_chr, all_start, all_end, all_rd, outlier_scores,label,index)
        #write_dataset_file( all_start,seg_rd,label)

        #lower = boxplot(outlier_scores)
        print("upper=",upper)

        CNV_chr, CNV_start, CNV_end, CNV_rd, CNV_type = combining_cnv(all_chr, all_start, all_end, all_rd, outlier_scores, upper,
                                                                  mode)
        cn = calculating_copy_number(mode, CNV_rd, CNV_type)
        write_cnv_file(CNV_chr, CNV_start, CNV_end, CNV_type, cn, outfile)
        sensitivity, precision, F1_score = calculate_scores(groudtruth, result_file, score_result_file)
        print("F1_score =", F1_score)
        print("precision =", precision)
        print("sensitivity =", sensitivity)
        end = datetime.datetime.now()
        print("running time: " + str((end - start).seconds) + " seconds")
        messagebox.showinfo("Completion", "Analysis is successfully completed!")
    except Exception as e:
        error_message = f"An error occurred: {str(e)}"
        print(error_message)
        messagebox.showerror("Error", error_message)

def browse_bam_file():
    bam_file_path = filedialog.askopenfilename(filetypes=[("BAM Files", "*.bam")])
    bam_file_entry.delete(0, tk.END)
    bam_file_entry.insert(0, bam_file_path)

def browse_ref_path():
    #ref_path = filedialog.askdirectory()
    ref_path =filedialog.askopenfilename(filetypes=[("BAM Files", "*.fa")])
    ref_path_entry.delete(0, tk.END)
    ref_path_entry.insert(0, ref_path)

def browse_output_file():
    output_file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text Files", "*.txt")])
    output_file_entry.delete(0, tk.END)
    output_file_entry.insert(0, output_file_path)

def browse_p_value_file():
    p_value_file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text Files", "*.txt")])
    p_value_file_entry.delete(0, tk.END)
    p_value_file_entry.insert(0, p_value_file_path)


    
def browse_groudtruth():
    groudtruth_path = filedialog.askopenfilename(filetypes=[("Groundtruth Files", "*.gt")])
    groudtruth_entry.delete(0, tk.END)
    groudtruth_entry.insert(0, groudtruth_path)

def browse_result_file():
    result_file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text Files", "*.txt")])
    result_file_entry.delete(0, tk.END)
    result_file_entry.insert(0, result_file_path)

def browse_score_result_file():
    score_result_file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text Files", "*.txt")])
    score_result_file_entry.delete(0, tk.END)
    score_result_file_entry.insert(0, score_result_file_path)

# Function to be executed when the "Run" button is clicked
def on_run_button_click():
    try:
        bam = bam_file_entry.get()
        ref_path = ref_path_entry.get()
        binSize = int(bin_size_entry.get())
        outfile = output_file_entry.get()
        p_value_file = p_value_file_entry.get()
        col = int(col_entry.get())
        k_value = int(k_entry.get())
        bandwidth = float(bandwidth_entry.get())  # Retrieve bandwidth from the entry widget
        groudtruth = groudtruth_entry.get()
        result_file = result_file_entry.get()
        score_result_file = score_result_file_entry.get()
        selected_method = segmentation_method.get()

        # Call the run_analysis function with all the parameters
        run_analysis(bam, ref_path, binSize, outfile, p_value_file, col, k_value, groudtruth, result_file, score_result_file, selected_method, bandwidth)

        print("Run button clicked with parameters:")
        print(f"BAM File: {bam}")
        print(f"Reference Path: {ref_path}")
        print(f"Bin Size: {binSize}")
        print(f"Output File: {outfile}")
        print(f"P-Value File: {p_value_file}")
        print(f"Column: {col}")
        print(f"k_value: {k_value}")
        print(f"bandwidth: {bandwidth}")
        print(f"groudtruth_path: {groudtruth}")
        print(f"result_file: {result_file}")
        print(f"score_result_file: {score_result_file}")
        print(f"Segmentation Method: {selected_method}")
    except Exception as e:
        error_message = f"An error occurred: {str(e)}"
        print(error_message)
        messagebox.showerror("Error", error_message)
# [I'm omitting them here to save space, but they should stay in your code]

# Create the main window with a modern theme
root = tk.Tk()
root.title("EDRAD-CNV Analysis for Real Data")
root.geometry("1000x800")  # More reasonable window size
root.minsize(900, 700)     # Set minimum window size

# Configure style
style = ttk.Style()
style.theme_use('clam')  # Modern theme

# Configure colors
bg_color = '#f0f0f0'
frame_color = '#ffffff'
accent_color = '#4b6cb7'
section_title_bg = '#3a7ca5'
section_title_fg = 'white'
text_color = '#333333'

style.configure('TFrame', background=frame_color)
style.configure('TLabel', background=frame_color, foreground=text_color, font=('Helvetica', 10))
style.configure('TButton', background=accent_color, foreground='white', font=('Helvetica', 10), borderwidth=1)
style.configure('TEntry', font=('Helvetica', 10), padding=5)
style.configure('TCombobox', font=('Helvetica', 10), padding=5)
style.configure('Title.TLabel', font=('Helvetica', 14, 'bold'), foreground=accent_color)
style.configure('SectionTitle.TLabel', font=('Helvetica', 10, 'bold'), 
               background=section_title_bg, foreground=section_title_fg, padding=5)
style.configure('ParamLabel.TLabel', width=18, anchor='e')  # Fixed width for parameter labels

# Main container with padding
main_frame = ttk.Frame(root, padding="20")
main_frame.pack(fill=tk.BOTH, expand=True)

# Title frame
title_frame = ttk.Frame(main_frame)
title_frame.pack(fill=tk.X, pady=(0, 20))
ttk.Label(title_frame, text="EDRAD-CNV Analysis Tool", style='Title.TLabel').pack()

# Create notebook (tabbed interface)
notebook = ttk.Notebook(main_frame)
notebook.pack(fill=tk.BOTH, expand=True)

# Input Parameters Tab
input_tab = ttk.Frame(notebook, padding=10)
notebook.add(input_tab, text="Input Parameters")

# File Input Section
file_section = ttk.Frame(input_tab)
file_section.pack(fill=tk.X, pady=5)

# Section title with colored background
file_title = ttk.Label(file_section, text="File Inputs", style='SectionTitle.TLabel')
file_title.pack(fill=tk.X, pady=(0, 5))

# File inputs container
file_container = ttk.Frame(file_section)
file_container.pack(fill=tk.X, padx=5, pady=5)

# BAM File
bam_frame = ttk.Frame(file_container)
bam_frame.pack(fill=tk.X, pady=2)
ttk.Label(bam_frame, text="BAM File:", style='ParamLabel.TLabel').pack(side=tk.LEFT, padx=5)
bam_file_entry = ttk.Entry(bam_frame)
bam_file_entry.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=5)
ttk.Button(bam_frame, text="Browse", command=browse_bam_file).pack(side=tk.LEFT, padx=5)

# Reference FASTA File
ref_frame = ttk.Frame(file_container)
ref_frame.pack(fill=tk.X, pady=2)
ttk.Label(ref_frame, text="Reference File:", style='ParamLabel.TLabel').pack(side=tk.LEFT, padx=5)
ref_path_entry = ttk.Entry(ref_frame)
ref_path_entry.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=5)
ttk.Button(ref_frame, text="Browse", command=browse_ref_path).pack(side=tk.LEFT, padx=5)

# Groundtruth File
gt_frame = ttk.Frame(file_container)
gt_frame.pack(fill=tk.X, pady=2)
ttk.Label(gt_frame, text="Groundtruth File:", style='ParamLabel.TLabel').pack(side=tk.LEFT, padx=5)
groudtruth_entry = ttk.Entry(gt_frame)
groudtruth_entry.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=5)
ttk.Button(gt_frame, text="Browse", command=browse_groudtruth).pack(side=tk.LEFT, padx=5)

# Output Section
output_section = ttk.Frame(input_tab)
output_section.pack(fill=tk.X, pady=5)

# Section title with colored background
output_title = ttk.Label(output_section, text="Output Files", style='SectionTitle.TLabel')
output_title.pack(fill=tk.X, pady=(0, 5))

# Output files container
output_container = ttk.Frame(output_section)
output_container.pack(fill=tk.X, padx=5, pady=5)

# CNV Output File
cnv_frame = ttk.Frame(output_container)
cnv_frame.pack(fill=tk.X, pady=2)
ttk.Label(cnv_frame, text="CNV Output File:", style='ParamLabel.TLabel').pack(side=tk.LEFT, padx=5)
output_file_entry = ttk.Entry(cnv_frame)
output_file_entry.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=5)
ttk.Button(cnv_frame, text="Browse", command=browse_output_file).pack(side=tk.LEFT, padx=5)

# P-value Output File
pval_frame = ttk.Frame(output_container)
pval_frame.pack(fill=tk.X, pady=2)
ttk.Label(pval_frame, text="P-value Output File:", style='ParamLabel.TLabel').pack(side=tk.LEFT, padx=5)
p_value_file_entry = ttk.Entry(pval_frame)
p_value_file_entry.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=5)
ttk.Button(pval_frame, text="Browse", command=browse_p_value_file).pack(side=tk.LEFT, padx=5)

# Result File
res_frame = ttk.Frame(output_container)
res_frame.pack(fill=tk.X, pady=2)
ttk.Label(res_frame, text="Result File:", style='ParamLabel.TLabel').pack(side=tk.LEFT, padx=5)
result_file_entry = ttk.Entry(res_frame)
result_file_entry.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=5)
ttk.Button(res_frame, text="Browse", command=browse_result_file).pack(side=tk.LEFT, padx=5)

# Score Result File
score_frame = ttk.Frame(output_container)
score_frame.pack(fill=tk.X, pady=2)
ttk.Label(score_frame, text="Score Result File:", style='ParamLabel.TLabel').pack(side=tk.LEFT, padx=5)
score_result_file_entry = ttk.Entry(score_frame)
score_result_file_entry.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=5)
ttk.Button(score_frame, text="Browse", command=browse_score_result_file).pack(side=tk.LEFT, padx=5)

# Parameters Section
param_section = ttk.Frame(input_tab)
param_section.pack(fill=tk.X, pady=5)

# Section title with colored background
param_title = ttk.Label(param_section, text="Analysis Parameters", style='SectionTitle.TLabel')
param_title.pack(fill=tk.X, pady=(0, 5))

# Parameters container - now vertical layout
param_container = ttk.Frame(param_section)
param_container.pack(fill=tk.X, padx=5, pady=5)

# Bin Size
bin_frame = ttk.Frame(param_container)
bin_frame.pack(fill=tk.X, pady=2)
ttk.Label(bin_frame, text="Bin Size (bp):", style='ParamLabel.TLabel').pack(side=tk.LEFT, padx=5)
bin_size_entry = ttk.Entry(bin_frame, width=20)
bin_size_entry.pack(side=tk.LEFT, padx=5)

# Column
col_frame = ttk.Frame(param_container)
col_frame.pack(fill=tk.X, pady=2)
ttk.Label(col_frame, text="Column:", style='ParamLabel.TLabel').pack(side=tk.LEFT, padx=5)
col_entry = ttk.Entry(col_frame, width=20)
col_entry.pack(side=tk.LEFT, padx=5)

# k Value
k_frame = ttk.Frame(param_container)
k_frame.pack(fill=tk.X, pady=2)
ttk.Label(k_frame, text="k Value:", style='ParamLabel.TLabel').pack(side=tk.LEFT, padx=5)
k_entry = ttk.Entry(k_frame, width=20)
k_entry.pack(side=tk.LEFT, padx=5)

# Bandwidth
bw_frame = ttk.Frame(param_container)
bw_frame.pack(fill=tk.X, pady=2)
ttk.Label(bw_frame, text="Bandwidth:", style='ParamLabel.TLabel').pack(side=tk.LEFT, padx=5)
bandwidth_entry = ttk.Entry(bw_frame, width=20)
bandwidth_entry.pack(side=tk.LEFT, padx=5)

# Segmentation Method
seg_frame = ttk.Frame(param_container)
seg_frame.pack(fill=tk.X, pady=2)
ttk.Label(seg_frame, text="Segment Method:", style='ParamLabel.TLabel').pack(side=tk.LEFT, padx=5)
segmentation_method = tk.StringVar(value="R")
method_menu = ttk.Combobox(seg_frame, textvariable=segmentation_method, 
                          values=["Python", "R"], state="readonly", width=18)
method_menu.pack(side=tk.LEFT, padx=5)

# Status/Log Tab (empty for now, can be implemented later)
log_tab = ttk.Frame(notebook)
notebook.add(log_tab, text="Execution Log")

# Run Button with better styling
run_frame = ttk.Frame(main_frame)
run_frame.pack(fill=tk.X, pady=10)
run_button = ttk.Button(run_frame, text="Run Analysis", command=on_run_button_click, style='TButton')
run_button.pack(pady=10, ipadx=20, ipady=5)

# Start the Tkinter event loop
root.mainloop()