import collections

import numpy as np

from proteomics import sequence


def getDEPeptides(prot_seq, max_expected_length,
                  min_expected_length=4, return_single_tryp=True):
    ''' For a given protein sequence, what are the expected LysC and Tryp peptides

    return_single_tryp = return LysC peptides with only a single Trypsin peptide
    '''

    lpep2tpep = collections.defaultdict(list)
    tpep2lpep = collections.defaultdict(list)
    tpeps = collections.defaultdict(list)

    for lysc_silica_pep in sequence.iteratePeptides(
            prot_seq, method='lysC', missed_cleavages=0,
            min_length=min_expected_length, max_length=float("Inf"), output_last=False,
            set_iso_to_leucine=True):

        if not "R" in lysc_silica_pep[0]:
            continue

        tryp_peps = sequence.iteratePeptides(
            lysc_silica_pep[0], method='trypsin', missed_cleavages=1,
            min_length=min_expected_length, max_length=max_expected_length, output_last=False)

        tryp_peps = list(tryp_peps)

        if len(tryp_peps) == 1 and not return_single_tryp:
            continue

        # if len == 0, loop not entered
        for tryp_silica_pep in tryp_peps:
            if tryp_silica_pep[0][0] == "M" and tryp_silica_pep[1] == 0:
                tryp_silica_pep_minus_meth = (
                    tryp_silica_pep[0][1:], tryp_silica_pep[1]+1,
                    tryp_silica_pep[2], tryp_silica_pep[3])

                lpep2tpep[lysc_silica_pep[0]].append(tryp_silica_pep_minus_meth)
                tpep2lpep[tryp_silica_pep_minus_meth[0]].append(lysc_silica_pep)
                tpeps[tryp_silica_pep_minus_meth[0]].append(tryp_silica_pep_minus_meth)

            lpep2tpep[lysc_silica_pep[0]].append(tryp_silica_pep)
            tpep2lpep[tryp_silica_pep[0]].append(lysc_silica_pep)
            tpeps[tryp_silica_pep[0]].append(tryp_silica_pep)

    return lpep2tpep, tpep2lpep, tpeps


def getDirectEvidence(df, prot2seq, outfile_name,
                      method="conservative",
                      EtOH_FT_threshold=1,
                      Kit_FT_threshold=1,
                      sequential_threshold=None,
                      threshold_method="both",
                      max_expected_length=100,
                      min_expected_length=6):

    #threshold_method = one of ["both", "either"]. both = use both thresholds, either = use either threshold
    
    # conservative method = ignore all LysC peptide where we see any Trpy peptide in EtOH control FT
    # mid-conservative method = ignore all Trpy peptide also seen in EtOH control FT
    # relaxed = ignore the EtOH control FT
    allowed_methods = ["conservative", "mid_conservative", "relaxed"]
    assert method in allowed_methods, "method must be one of: %s" ",".join(allowed_methods)
    
    
    outfile = open(outfile_name, "w")
    outfile.write("%s\n" % "\t".join(map(str, (
        "uniprot_id", "lysC_pep_ix", "lysC_pep_start", "lysC_pep_end",
        "binding_site_type", "site_start", "site_end", "binding_seq",
        "method", "EtOH_FT_threshold", "Kit_FT_threshold"))))

    for uniprot_id in set(df['master_protein']):
        #uniprot_id = 'A2RU67'
        single_prot_df = df[df['master_protein']==uniprot_id]
        if threshold_method == "both":
            if method in ["conservative", "mid_conservative"]:
                test_df = single_prot_df.loc[(single_prot_df["EtOH_FT"] >= EtOH_FT_threshold) &
                                             (single_prot_df["Tryp_FT"] == 0 ) &
                                             (single_prot_df["Kit_FT"] >= Kit_FT_threshold) &
                                             (single_prot_df["KitTryp_FT"] == 0),:]
            else:
                test_df = single_prot_df.loc[(single_prot_df["EtOH_FT"] >= EtOH_FT_threshold) &
                                             (single_prot_df["Kit_FT"] >= Kit_FT_threshold), :]
        elif threshold_method == "either":
            if method in ["conservative", "mid_conservative"]:
                test_df = single_prot_df.loc[((single_prot_df["EtOH_FT"] >= EtOH_FT_threshold) &
                                              (single_prot_df["Tryp_FT"] == 0)) |
                                             ((single_prot_df["Kit_FT"] >= Kit_FT_threshold) &
                                              (single_prot_df["KitTryp_FT"] == 0)), :]
            else:
                test_df = single_prot_df.loc[(single_prot_df["EtOH_FT"] >= EtOH_FT_threshold) |
                                             (single_prot_df["Kit_FT"] >= Kit_FT_threshold),:]
        else:
            raise ValueError("threshold method must be one of ['both', 'either']")

        if sequential_threshold:
            test_df = test_df.loc[test_df['Sequential_FT'] == 0, :]

        prot_seq = prot2seq[uniprot_id]
        
        lpep2tpep, tpep2lpep, tpeps = getDEPeptides(
            prot_seq, max_expected_length, min_expected_length)

        lpep_ignore = set()

        if method == "conservative":
            false_positive_df = single_prot_df[
                (single_prot_df["Tryp_FT"] > 0) | (single_prot_df["KitTryp_FT"] > 0)]
            for pep_seq in false_positive_df['Sequence']:
                pep_seq = pep_seq.replace("I", "L") 
                if pep_seq in tpep2lpep.keys():
                    lpep_ignore.update(tpep2lpep[pep_seq])

        lpep_hits = set()
        for pep_seq in test_df['Sequence']:
            pep_seq = pep_seq.replace("I", "L")
            if pep_seq in tpep2lpep :
                lpep_hits.update(set(tpep2lpep[pep_seq]))
                
        # remove lpep also identified from false positive Trypsin peptides
        # note, if method!="conservative", lpep_ignore will be an empty set
        lpep_hits = lpep_hits.difference(lpep_ignore)

        prot_hits = np.zeros(len(prot_seq), dtype=bool)

        for lpep in lpep_hits:
            prot_hits[lpep[1]:lpep[2]] = 1

        single_tryp_cover = set()
        for pep_seq, start, end in test_df[
            ['Sequence', 'peptide_start', 'peptide_end']].itertuples(index=False):
            pep_seq = pep_seq.replace("I", "L")
            if pep_seq in tpep2lpep:
                for lpep in tpep2lpep:
                    if lpep == pep_seq:
                        # tryp peptide covers the whole LysC peptide (should we remove these)?
                        # right now we're just marking these for removal later
                        single_tryp_cover.add(lpep)
                poss_t_peps = tpeps[pep_seq]

                if len(poss_t_peps) == 0:
                    raise ValueError("peptide sequence is not sequential - shouldn't arrive here")
                elif len(poss_t_peps) > 2:
                    print("too many possible locations - need to work out what to do here", uniprot_id)

                for t_pep in poss_t_peps:
                    t_pep_seq, tpep_start, tpep_end = t_pep[0:3]
                    poss_l_peps = tpep2lpep[t_pep_seq]

                    if len(poss_l_peps) == 0:
                        raise ValueError("this shouldn't happen!")

                    for l_pep in poss_l_peps:
                        prot_hits[tpep_start+l_pep[1]: tpep_end+l_pep[1]] = 0

            else:
                pass


        for ix, l_pep in enumerate(lpep_hits):
            pep_aa_hits = prot_hits[l_pep[1]: l_pep[2]]

            
                
            
            if l_pep[1] == 0 and pep_aa_hits[1] == 0:
                pep_aa_hits[0] = 0

    
            if sum(pep_aa_hits)== 0:
                if l_pep[0] in single_tryp_cover:
                    binding_site_type = "Single Tryp Pep Cover"
                else:
                    binding_site_type = "Multiple Tryp Pep Cover"

                site_start, site_end = map(int, l_pep[1:3])
                binding_site_seq = prot_seq[site_start:site_end+1]
                
                outfile.write("%s\n" % "\t".join(map(str, (
                    uniprot_id, ix, l_pep[1], l_pep[2], binding_site_type, site_start, site_end+1,
                    binding_site_seq, method, EtOH_FT_threshold, Kit_FT_threshold))))

                continue

            transitions = sum(pep_aa_hits[1:] ^ pep_aa_hits[:-1])
            if transitions == 1 or pep_aa_hits[0] == False and transitions == 2: # 2 transitions OK if starts with 0
                binding_site_type = "unique"
            elif transitions > 1:
                binding_site_type = "non-unique"
                site_start, site_end = map(int, l_pep[1:3])
                binding_site_seq = prot_seq[site_start:site_end+1]
                
                outfile.write("%s\n" % "\t".join(map(str, (
                    uniprot_id, ix, l_pep[1], l_pep[2], binding_site_type, site_start, site_end+1,
                    binding_site_seq, method, EtOH_FT_threshold, Kit_FT_threshold))))
                continue
                
            else:
                print(l_pep)
                print(pep_aa_hits)
                print(test_df)
                raise ValueError("unexpected value %s: %s" % (transitions, uniprot_id))

            binding_site_positions = np.array(range(l_pep[1], l_pep[2]))
            binding_site_positions = binding_site_positions[pep_aa_hits]


            binding_sites = []
            current_site = [binding_site_positions[0],binding_site_positions[0]]
            for pos in binding_site_positions[1:]:
                if pos != current_site[1] + 1:
                    binding_sites.append(current_site)
                    current_site = [pos, pos]
                current_site[1] = pos
            binding_sites.append(current_site)

            for binding_site in binding_sites:
                site_start, site_end = binding_site
                binding_site_seq = prot_seq[site_start:site_end+1]
                outfile.write("%s\n" % "\t".join(map(str, (
                    uniprot_id, ix, l_pep[1], l_pep[2], binding_site_type, site_start, site_end+1,
                    binding_site_seq, method, EtOH_FT_threshold, Kit_FT_threshold))))

    outfile.close()


def makeProtterInput(df, uniprot_ids, outfile=None):
    ''' make output which can be loaded into protter

    df = direct_evidence pandas dataframe. This should be filtered
       beforehand to exclude unwanted sites. Must contain columns:
       - 'uniprot_id'
       - 'site_start'
       - 'site_end'
    uniprot_id = uniprot_ids of interest
    outfile = filename for outfile, if None, prints out
    '''
    if outfile:
        out = open(outfile, "w")
    
    for uniprot_id in uniprot_ids:#.intersection(consistent_TMT):
        tmp_df = df[df['uniprot_id']==uniprot_id]
        
        if outfile:
            out.write("%s\t" % uniprot_id)
        else:
            print(uniprot_id, end="\t")
        
        for start, end in tmp_df[['site_start', 'site_end']].itertuples(index=False):
            start, end = map(int, (start+1, end))
            if outfile:
                out.write("%s-%s," % (start,end))
            else:
                print("%s-%s," % (start,end), end="")
            
        if outfile:
            out.write("\n")
        else:
            print("\n", end="")
    
    if outfile:
        out = open(outfile, "w")
