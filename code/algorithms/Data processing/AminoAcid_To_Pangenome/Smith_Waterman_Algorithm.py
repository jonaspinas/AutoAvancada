
import itertools
import numpy as np

#Match score and gap cost can be changed and reworked to better fit the model
#For this algorithm match score and gap cost were kept at 1 and -2 (most common scoring system across multiple searches)
#traceback and smith waterman are mainly visualization tools, the algorithm is contained in the matrix function
#For a more in-depth explanation of the algorithm follow the data treatment section of the delivery
def matrix(a, b, match_score=1, gap_cost=2):
    H = np.zeros((len(a) + 1, len(b) + 1), int)

    for i, j in itertools.product(range(1, H.shape[0]), range(1, H.shape[1])):
        match = H[i - 1, j - 1] + (match_score if a[i - 1] == b[j - 1] else - match_score)
        delete = H[i - 1, j] - gap_cost
        insert = H[i, j - 1] - gap_cost
        H[i, j] = max(match, delete, insert, 0)
    #The score of the algorithm is the max of the matrix
    Score = np.max(H)
    Average_len = (len(a)+len(b))/2
    #To get a gene similariy in percentage, the score is normalized by the average gene length 
    Normalized_Score = Score/(Average_len * match_score)
    return H, Score, Normalized_Score

def traceback(H, b, b_='', old_i=0):
    # flip H to get index of last occurrence of H.max() with np.argmax()
    H_flip = np.flip(np.flip(H, 0), 1)
    i_, j_ = np.unravel_index(H_flip.argmax(), H_flip.shape)
    i, j = np.subtract(H.shape, (i_ + 1, j_ + 1))  # (i, j) are last indexes of H.max()
    if H[i, j] == 0:
        return b_, j
    b_ = b[j - 1] + '-' + b_ if old_i - i > 1 else b[j - 1] + b_
    return traceback(H[0:i, 0:j], b, b_, i)

def smith_waterman(a, b, match_score=1, gap_cost=2):
    a, b = a.upper(), b.upper()
    H, score ,nscore = matrix(a, b, match_score, gap_cost)
    b_, pos = traceback(H, b)
    return pos, pos + len(b_), score, nscore

#Uncomment the below code for an example of the working code
'''
a, b = 'hipopotamo', 'hipoppotamo'
start, end, score, nscore = smith_waterman(a, b)
print(a[start:end], score, nscore)     # GTTGAC
'''