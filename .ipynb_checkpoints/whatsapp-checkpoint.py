import operator
from matplotlib import pyplot as plt
    
def variations_of_word(corpus, variation_list):
    """
    Finds words in corpus that start with any of the words given in the variation_list
    
    Parameters
    ----------
    corpus : dict
    variation_list : list of str
    
    Returns
    -------
    subset_dict: dict
    """
    subset_dict = {}
    for var in variation_list:
        for word in corpus.keys():
            if word.startswith(var):
                subset_dict[word] = corpus[word]
    return subset_dict

def sort_dict_by_values(adict, method='descending'):
    """
    Sort adict by values
    """

    if method not in ['descending', 'ascending']:
        raise ValueError('method must be descending or ascending')
    sorted_x = sorted(adict.items(), key=operator.itemgetter(1))
    if method == 'ascending':
        return sorted_x
    else:
        return list(reversed(sorted_x))
    
def count_user_messages(lines, names_dict):
    
    # Set all counters of names in dict to 0
    for name in names_dict.keys():
        names_dict[name] = 0
    
    for line in lines:
        for name in names_dict.keys():
            if name in line:
                names_dict[name] += 1
    return names_dict
                
def get_all_lines(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()
    return lines

def plot_messages_per_user(names_dict, ax=None):
    if ax is None:
        f, ax = plt.subplots(figsize=(7, 5))
    
    sorted_names = sort_dict_by_values(names_dict)
    
    names = [x[0] for x in sorted_names]
    values = [x[1] for x in sorted_names]
    
    ax.bar(range(len(names)), values, align='center')
    ax.xticks(range(len(names)), names)
    ax.ylabel('Messages')
    
    return ax
    