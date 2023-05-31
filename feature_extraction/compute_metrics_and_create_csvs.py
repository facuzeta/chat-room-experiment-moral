import pandas as pd
import json
import sys
from tqdm import tqdm


def compute_df_s1(data):
    df_s1 = pd.DataFrame([r for d in data for r in d['stage_1']])
    dic_participant_id__group_id = {
        r['participant_id']: d['group_id'] for d in data for r in d['stage_1']}

    res_s1 = []
    for p_id in tqdm(df_s1.participant_id.unique()):
        r = {}
        r['participant_id'] = p_id
        r['group_id'] = dic_participant_id__group_id[p_id]
        df_subj = df_s1[df_s1.participant_id == p_id]
        for i, qid in enumerate(range(9, 15)):
            if len(df_subj[df_subj.question_id == qid]) == 0:
                continue
            r[f'q{i+1}answer'] = df_subj[df_subj.question_id ==
                                         qid].answer_value.values[0]
            r[f'q{i+1}confidence'] = df_subj[df_subj.question_id ==
                                             qid].confidence_value.values[0]
        res_s1.append(r)
    return pd.DataFrame(res_s1)


def compute_df_s3(data):
    
    df_s3 = pd.DataFrame([r for d in data for r in d['stage_3']])
    dic_participant_id__group_id = {
        r['participant_id']: d['group_id'] for d in data for r in d['stage_3']}

    res_s3 = []
    for p_id in tqdm(df_s3.participant_id.unique()):
        r = {}
        r['participant_id'] = p_id
        r['group_id'] = dic_participant_id__group_id[p_id]
        df_subj = df_s3[df_s3.participant_id == p_id]
        for i, qid in enumerate(range(9, 15)):
            if len(df_subj[df_subj.question_id == qid]) == 0:
                continue
            r[f'q{i+1}answer'] = df_subj[df_subj.question_id ==
                                         qid].answer_value.values[0]
            r[f'q{i+1}confidence'] = df_subj[df_subj.question_id ==
                                             qid].confidence_value.values[0]
        res_s3.append(r)

    df_res = pd.DataFrame(res_s3)   
    return df_res


def compute_metrics_chats(data):
    """
    Compute:
     - number of words used by participant in question
     - number of interventions by participant in question
     - Mean word length of words sent by this participant in question
    """

    res_features_by_participant = {}
    for d in tqdm(data):
        for i in range(1, 5):
            stage = f'stage_2_{i}'
            df_chat = pd.DataFrame(d['stage_2'][stage]['chat'])
            if len(df_chat) == 0: continue

            # reindex 9->1, 10->2, ...
            question_id_reindexed = d['stage_2'][stage]['question_id']-8

            for p_id in df_chat.participant_id.unique():
                if p_id not in res_features_by_participant:
                    res_features_by_participant[p_id] = {'participant_id': p_id}
                df_chat_one_subject = df_chat[df_chat.participant_id == p_id]

                res_features_by_participant[p_id][f'q{question_id_reindexed}wordsTotal'] = compute_number_of_words(df_chat_one_subject)
                res_features_by_participant[p_id][f'q{question_id_reindexed}interventionsTotal'] = compute_number_of_interventions(df_chat_one_subject)
                res_features_by_participant[p_id][f'q{question_id_reindexed}wordLengthMean'] = compute_mean_word_length(df_chat_one_subject)
    
    df_features = pd.DataFrame(res_features_by_participant.values())
    return df_features


def compute_number_of_words(df_chat_one_subject):
    """
    Compute number of words used by participant in question
    """
    n_words = sum(
        [1 + e.count(" ") for e in df_chat_one_subject.text]
    )
    return n_words

def compute_number_of_interventions(df_chat_one_subject):
    """
    Compute number of interventions by participant in question
    """
    n_interventions = len(df_chat_one_subject)
    return n_interventions

def compute_mean_word_length(df_chat_one_subject):
    """
    Compute Mean word length of words sent by this participant in question
    """
    words_length_series = [len(w) for l in df_chat_one_subject.text.tolist() for w in l.split()]
    mean_word_length = sum(words_length_series) / len(words_length_series)
    return mean_word_length





def compute_and_create_csvs(data):
    """
    compute_and_create_csvs
    """

    clean_data = [d for d in data if len(d['stage_2'].keys()) == 4]

    df_s1 = compute_df_s1(clean_data)
    df_s3 = compute_df_s3(clean_data)
    df_features = compute_metrics_chats(clean_data)

    df_s1.merge(df_features, how='outer').to_csv('s1.csv', index=False)
    df_s3.merge(df_features, how='outer').to_csv('s3.csv', index=False)


if __name__ == "__main__":

    # use first argument as input file
    fn = sys.argv[1]
    print("Using file", fn)

    # read data
    with open(fn, "r") as f:
        data = json.load(f)
    print("Read", len(data), "groups")

    # create csvs
    compute_and_create_csvs(data)
