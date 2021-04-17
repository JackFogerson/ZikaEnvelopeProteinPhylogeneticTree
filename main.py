from Bio import Phylo
from io import StringIO
import numpy

#compares sequences and adds to matrix
def sequenceCompare():
    #original sequence (Brazilian Zika Strain Envelope Protein)
    sequenceBrazil =        'GGGACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAGTCCAAG'

    #25 Basic Data Sets (Regional Zika Strain Envelope Protein)
    sequenceVietnam =       'GGGACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAGGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAGTCCAAA'
    sequencePuertoRico =    'GGGACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACTGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAGTCCAAA'
    sequenceGuatemala =     'GGGACTCCACACTGGAACAACAAAGAAGCACTGGTGGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAGTCCAAA'
    sequenceGuadeloupe =    'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAGCCCA  '
    sequenceUSA =           'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAGCCCA  '
    sequenceMexico =        'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAGCCCA  '
    sequenceMartinique =    'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAGCCCA  '
    sequenceYucatan =       'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAG      '
    sequenceHaiti =         'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAG      '
    sequenceFiji =          'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAG      '
    sequenceFrenchGuiana =  'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAG      '
    sequenceAngola =        'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAG      '
    sequenceCuba =          'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAG      '
    sequenceHonduras =      'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACNNGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAG      '
    sequenceRioDeJaneiro =  'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAANGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAG      '
    sequenceColombia =      'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAG      '
    sequenceBoracay =       'GGGACTCCACATTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTCAGATTGAAGGGCGTGTCATACTCCTTGTGTACTGCAGCATTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAGTCCAAA'
    sequenceRioDeJaneiro2 = 'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACYGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGKTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAG      '
    sequenceJamaica =       'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAKGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGKACCTTGCAAGGTTCCAG      '
    sequenceDominicanRep =  'GGAACTCCACACTGGAAYAACAAAGAAGCACTGGTAGAGTTCAGGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAG      '
    sequenceCuiaba =        '                                                                                                                                       ATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGTTTTGTCCCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGTTTCACTGCAGGGCACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCTTGCAAGGTTCCAGCCCAA '
    sequenceIndia =         '      CCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCCGGAGCTCTGGAGGCTGAGATGGATGGTACAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGAT                                                                                   '
    sequenceBeloHorizonte = 'GGAACTCCACACTGGAACAACAAAGAAGCACTGGTAGAGTTCAAGGACGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGTTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACTTAGATTGAAGGGCGTGTCATACTCCTTGTGTACTGCAGCGATCACATTCACCAAGATCC                                                                                 '
    sequenceIndia2 =        'GGAACTCCACATTGGAACAACAAAGAAGCATTGGTAGAGTTCAAGGATGCACATGCCAAAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCCGTTCACACGGCTCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGTTGTCTTCTGGCCACTTGAAATGTCGTCTGAAAATGGATAAACTTAGATTGAAGGGTGTGTCATACTCCTTGTGCACTGCGGCGTTCACATTCACCAAGATCCCGGCTGAAACGCTGCATGGGACAGTCACAGTGGGGGTACAGTA                                      '
    sequenceMexico2 =       '                                                      GCCAGAAGGCAAACTGTCGTGGTTCTAGGGAGTCAAGAAGGAGCAGCTCACACGGCCCTTGCTGGAGCTCTGGAGGCTGAGATGGATGGTGCAAAGGGAAGGCTGTCCTCTGGCCACTTGAAATGTCGCCTGAAAATGGATAAACCTAGATTGAAGGGCGTGTCATACTCCTTGTGTACCGCAGCGTTCACATTCACCAAGATCCCGGCTGAAACACTGCACGGGACAGTCACAGTGGAGGTACAGTACGCAGGGACAGATGGACCT                   '

    #25 Related Data Sets (Other Viral Envelope Proteins)
    sequenceDNUSBat =       '----------------------AAGAAGCAGAGAA---GTT-------------GC----GGAGGGAGTGTCTTTCTTCTGTGGTTTTGTGTTGGCTTACTTGGTTCTGTCATC--AGGAGTCATGGGGGGTGGACTATATCA---CCATCCAAGGCTGAGGTTGGCCAACGCTAAAGGTTGGAAATGGATCGTGCATCTTTGCCTCCCTGGACATTGGGAGTCCATGTGAAAGGACCATTTCTTA-------------------------TGAATGTGTTACTCTGAGTGAAAAAG----------------AAGAACCATTTGACGTTGAT-------'
    sequenceHIVNetherlands ='---------------------------------------------------------------------GTAGTA---ATTAGATCT---------GACAATTTCACGGACAATACTAAAACCATAATAGTACAGCTGAACAAATCTGTAGAAATTAATTGCACAAGACCCA-ACAACACCACAAGAAAAAGTATACATATAGGACCAGG------G--AGAGCATTTTATACACCAGGAGAAATAAT-------------------AGGAGATATAAGACAAGC--------------------------------ACATTGT----------------'
    sequenceSepticemia =    '-----------GTGG--CAGGCCATCATCCCT------------GGGGACTCACAATGGCATGCACAGTGACATTTTGCGGGACAGATT----G-----GATCAAGACTGACCTGGGGGATCTGATCCAGGTGACAGGACCGGAGGGCAGAGGAAACTGAC-TCCAAAC-----AAGTGTGT---CAACACCGATGT---CCAGATGAGGGGGGCAACAGACGACTTCTCTTTCTCAACCATCTCAT-CACCAACATGGCTCAAAAACCGAGTGCCTAGATGCCCATAGTGATATCACCGCTTCTGGGAAAGTCTCCTCATTTCT---------------'
    sequenceXMRV =          '----GAGACCACTGGA--------------CAGGCATACTG----GAAGCCATCATCATCATGGGACCT----------------------------------------------------------------------------------------------------------AATTTCCCTTAAGCGAGGAAACACTCCTAACAGGGCCCCTGTTATGATTCCTCGGTC--------------TCCAGTGGCGTCCAGGGTGCCACACCG--------------------------------GGGGGTCGAT-------GCAACCCCCTAGTCTTA-'
    sequenceHepCEgypt =     '----------------------------------------------------ATCACCTACCGCAATGCCTCGGGCATCTACCATGTCACCAATGACTGC----------CCGAACTCAAGCATAG-TATATGAGGCCGACCATCACATCCTACA---TCTCCCAGGCTGCGTGCCTTGT-GTGAGAATGGGGAATCAGTCACGTTGTTGGGTGGCCCTTACTCCTACCGTCGCAG-----CGCCATACATCGGCGCGCCGCTTGAGTCCCTGCGGAGTCATGTGGATCTGATGGTGGG-----------------------------------------'
    sequenceHepCFrance =    'ATGATGATGAACTGGTCGCCCACAACAGCCATAGTGGTGTCGCAGTTACCGGATCCCACAAGCCATCGTGGACATGG--TGGCGGGGGCCCACTGGGGAGTCCTGGCGGGCCTCGCCTACTATTCCAGAACTGGGCTAAGGTCT-TAATTGTGATGCTGCTCTTTGCCGGCGTTGAC--------------GGGGAGACCACCGTGGTGGGGGGAAGCCAGATCCATCCGTAGTGTAACATCCCTCTTTAGTTTTGGGGCATCTCAGAAAATCCAGCTCATAAACAC-----------------CAACGGCAGCTGGCACATCAACAGGACTGC------'
    sequenceHIVJapan =      'TGTACAAGACCCGGCAACCATACAAGAAAAGTAAT-------------ACATATAGGACCAGGGAGAGTATATTATACAACAGGAGAA---------------------------------------------------------------------------------------------------------ATAATAGGAGATATAAGACGA------GCACATTGT---------------------------------------------------------------------------------------------------------------'
    sequenceHIVItaly =      'CTAAAATCATAATAGTACAGCTGAATAACACTGTAAATATTACTT------------------------GTATAAGACCCAACAACAATACAAGAAAAAGTACCCATGGGACCAGGGAGAGCATTTTAT---------------GCAACAGGAGACATAAT-AGGAGATATAAGAAAAGCACAT---------TGTAACATTAGTAGAATAGACTGGAATAACACTTTAAAACAGATAGC-----------------------TGAAAAATTAAGA-GAACAATTT------------------------------------------------------'
    sequenceDengue =        '----CTCTTCACTGGACATCTCAA---------GTGCAGGCTGAGAATGGACAAGCTACAGCTCAAAGGAATGTCATACTCTACACAGGAAAGTTTAAAGTTGTGAAGGAAATAGCAGAAACACAACATGGAACGATAGAGAGTGCAATATGAGGGGACGGCTCTCCATGCAAAATCCCTTTCGAGATAATGGATT----T-GGAAAAAAGATATGTCTTAGGCCACC----------TGATCACAGTCAACCCAATTGTAACAGAAAAAGA---TAGCCCAGTCAACATAGAAGCAGAACATTCGGAGACAGTTAC-----------------------'
    sequenceDengue2 =       '---------------------ATAGGAGTCATCATCACATGGATGGAATGATTCACGTAGCACCTCACTGTCTGTGTC---ACTCGTATTGGTGGGAGTCGTAACACTGTACTTG------------------GGAGTCATGGTGCAGGCTGATAGT------------GGTTGCGTTGTGAGCTGGAAAAACAAAAACTGAAATGTGGCAGCGGG-----ATTTTTATCACAGATAACGTGCACACA-----------------------TGGACAGAACAATATAAGTTCCA------A--------CCAGAATCCCCTTCAAAACTAGCT-------'
    sequenceEncephalitis =  '--------------------------TACAGTGA-----TAGGAGAACACGCCTGGGACTTTGGTTCCACTGGGGGTTTCTTGGCTTGTGGGCAAAGCGCTGCACACAGTTCTCGGTGGCGCCTTCAACAGCATCT---TTGGGGGAGTTG-GGTTTTTTCCCAAGCTCCTGTTGGGTGTAGCCTTGGCTTGGTTGGGCCT--------------------GAACATGAGGAACCCCACCATGTCCATGAGCTTCCTCCTGGCTGA--------------------------------------------------------------------------'
    sequenceOrf =           'CGAACTTCCACCTCGACC-ACTCCGGGGGCGGCGTATTCTTCTCGGACTCGC----CGGAGCGCTTCCTAGGCTTCTACCACCGCATCGAGAACGCCAAGAACAGCATCGACCTCTCGCTGCTCTCGATGGTGCCGGTGATCAAGC------ACGCCAGCGC-----CGTGGAGTACTGGCCGCGGATCATAGACGCACTGCTGCGCGCGG--CCATCGACCGCGGCGTGCGCGTGCGCGTGATCATTACCGAGGACCCGCTGTCGGTCTCGGCCGCGCGCAGCCTCGACGACTTTGGCGTCGGCAGC-----GTGGACATGTCCGTGCGCAAGTTCGTG'
    sequenceOrf2 =          '--AACTTCCACCTCGACC-ACTCCGGGGGCGGCGTATTCTTCTCGGACTCGCCTCCTAGGCTTCTACAGGACCTCGTTGCACCGCATCGAGAACGCCAAGAACAGCATCGACCTCTCGCTGCTCTCGATGGTGCCGGTGATCAAGCA------CGCCAGCGC-----CGTGGAGTACTGGCCGCGGATCATAGACGCACTGCTGCGCGCGG--CCATCGACCGCGGCGTGCGCGTGCGCGTGATCATTACCGAGGACCCGCTGTCGGTCTCGGCCGCGCGCAGCCTCGACGACTTTGGCGTCGGCAGCG-----TGGACATGTCCGTGCGCAAGTTCGTG'
    sequenceHerpes =        'AATGCTTCACAATCAGACATATGAAAAGTACGGAAACGTGTCGTCTTCGAAACTA----CCGGCGGACTAGTAGTGTTCTGGCAAGGTATCAAGCAAAAAT---------CTCTGGTGGAACTCGAACG--TTTGGCCAACCGTTCCAGTCTGAATCTTACTCATAACCAGTACAGATGGCCAATGTAACTCATTTATCTAATATGGATTCGGTAC---ACAATCTGGTCTACGCCCA-GCTGCAGTT---------CACCTACGACACGTTGCGCGGCTACATCAACCGGGCGTTGAC-----------------------------------------'
    sequenceSheeppox =      '------GCTCTTTTGT---------------------------ACAAAGAGCATTACATAATCCAGAAAAATATTC----------------TGTAAAATTTTCAACACCTCCTGATTTTTCTACCTTTTCCCATATAAGTAACTTATATGATAAACTGATATTTTATCTTTAAAATTTACATCTGAATTTTTAAAATCTTTTACT-------GTGTCAACTTTTTTATAAAATATATCATTGTCACTTTTTAATTCTGGAACTACATCTGAAATTTCGCGACCAACGATTGGTATAACATATAATGGGATATCTGCCATTTGATAATTAAAAATTAAAG'
    sequenceWhiteSpot =     'AACACTGGCCCTTCTTACACCATGGAAGATCTTGAAGGCTCCATGTCTATGGCTCGCATGGGTCTC----TTTTTGATCGTTGCTATCTCAATTGGTATCCTCGTCCTGGCCGTCATGAATG------------TATGGATGGGACCAAAGAAGGACAGCGATTCTGACACTGATAA-------------------------------------------------------GGACACCGATGATGATGACGACACTGCCAACGATAACGATGATGAGGACA-----AATATAAGAACAGGACCAGGGATATGATGCTTCTGGCTGGGTTTCCGCCCAAG'
    sequenceEColi =         '--------------------------------CTCACAGTTGACGGAAACGGCTACGTAACGTGAAAGTATCGGAAGTAACGCGTTTCACTTGCCGCCAGTTCACGGTCGACTTCATCCAGCAATTCGCTGCTGGAAGGTAGCTGCTC------------------ACCAAAGGCCATGTTGTCCGGCATTGGATCAAGCGCAAATGAC-------------------------------------------AAAATGATGTAAGCA--------------------------ACCAGGGTAAACAGCGCCAGACCGAA---------------------'
    sequenceNewCastle =     '------------------------------------------------------------------ACTACATCCGG---AGGGGGGAGACAGGGGCGCCTTATAGGCGCCATTATTGGCGGTGTGGCTCTTGGGGTTGCAACTGCCGCACAAATAACAGCGGCCGCAGCTCTGATACAA--------------GCCAAACAAAATGCTGC---CAACATCCTCCGACTTAAAGAGAGCATTGCCGCAAC------------------CA------------AT-------GAGGCTGTGCATGAGGTCACTGACGGATTATCCCAACTAGCG-------'
    sequenceHanta =         '---TTGTCAGTTTGGTGATTGGTGATATCATGAGTACACCTAGAGGGATGCAGTGCCCTGAA---CATGCTGGTTCATTTCGGAAG-AAGTGTGCATTTGCTACTACACCAGTATGTGAATATAGTGGGAATACAATCTCTGGCTATCAGAGGATGCTTGCAATCAGTCATTTAACATAACAGAACCACATAT----------AACTA---ATAACTCATTAGACTGGGTTGACCCAGATAGCTCATTAAAGGATCACATCAACTTAGTTGTTAATAGAGATGTCTCATTCCAAGACTTAT---------CTGAAAATCCCTGTCAAGTTGGTGTGGCTG'
    sequenceIslaVista =     '--------------------------------CTGTCAGTTCGGGGACCCAGGTGATATAATGAGTACCCCACAAGGCATGAGGTGTCCTGAACGAAATGTTCCCACAAGCCCTCATGCCAGTATGAAGGATACCATTTCTGGATATAAGAGAATGATAGCTACCAATCATTTAATACAGAACCACACATTACAGCAAATTCGCTTAAGGACCATGTAAATGTCCTAGTCAATAGAGATATTTCCTTTCAAGATCT----------------------------------------------------TTCTGAGACACCTTGTAAAGTTGATCTTTCAG'
    sequenceAnemiaVirus =   '-----------------------------------------------------------CAAGTTAATCTTAGTTCCTCT----------------------------AACTCCTCTGTACGTGTTGAGGAT---------GTAGCGAACACAAC---GGAATATTGGGGATTTAAATGGCTA--GAATGTAATCAAACAGAAAATTTTAAGACTATATTAGTACCTGAAAATGAAATGGTAAATATCAATGATAATGATACCTGG--------------------------------------------------------------------------'
    sequenceAnemiaVirus2 =  'GGGGTTCCTTCCCTGTAGACCATTTCAAAATTACTTCAGTTGAAGCATGCATATGGATAATA--ATACTGCTACATTATTAGAAGCTTATCATAGA-GAGATAACA--TTCATTTATAAGTC--TTCTTGCAC----AGATAGTGATCA------------TTGTCGAGA------GTATCAATGTAAAAAAGTTAATCTTAATTCCTCTGACTCCTCTAACTCTGTACGTGTTGAGGATGTAACGAACAGCGATATTGGGGATTTAAGGCTAGAATGTAATCAAACAGAAAATTTTAAGACTATATTAGTACCTGAAAATGAAATGGTAAATATCAAGG'
    sequenceSIV =           'CGTATTCTTAACTTTATTGGCCCCTGGGCATCCCAGTGCTTCAGGCCAGCCGGTGTACAATCACGGTAGGTATCTCCTCCTACCACTCCAGCCCCTGCAGCCCAGCTCAGCCTTTATGTACCTGGGCCCTCGACCTTGTGTCCATCACTAAGGACCAGCTCCTCTACCCCCCC-------TGCCAAAACCTGATCACCTATTCCAACTACCACAAGACCTACTCCCTGTATCTCTTCCCACAC------------------------TGGGTACAAAAGCCACTCCGCCGGGGGCTTTACTACTCAGCCTACTCTGATCCTTGCTCCCTGCAA-------'
    sequenceSiv2 =          '-------------------------GACAGATGGGGGTTAACAAAGAATGAAATACCAACAACAACAGCATCAACAAAACAAAAAGCCCAAAAGCAGAAGCAATAACAGCAAAAG---------TTATAAATGAA----------------------------------------------------------------------------------AGTGATCCTTGTATAAGCAACAATAA----------TTGTACAGGCTTAGAGCAGGAACCAATGGTA---------------------AGCTGTAGATTTAACATGACAGGGTTAA-------'
    sequenceHorseHerpes =   '------------ATGAACCTTAACGAT---GTGCC----CACCCTGCACACCATGATCACCCTGAACCTGTCTCTGGTGGAAAATGTAGATTTCCAGGTGATACAGCTTTACTCTCAAAAGGAAAAGAAGCTGTCTAATGTGTTTGACATAGAGACCATGTTTAGGGAGTACAACTACTACACCCAAAATCTAAAGGGGCTGAGA---AAGGATCTGGATGACTCCATCCACGATGGCAGGGACAGCTTTATCCAGTTCTTGGGGGACCTGGTGCAGGAC-------------------------------------------------------------'

    # add sequences to array and get length
    sequence_list = [sequenceBrazil, sequenceVietnam, sequencePuertoRico, sequenceGuatemala,sequenceGuadeloupe,
                     sequenceUSA, sequenceMexico, sequenceMartinique, sequenceYucatan, sequenceHaiti, sequenceFiji,
                     sequenceFrenchGuiana, sequenceAngola, sequenceCuba, sequenceHonduras, sequenceRioDeJaneiro,
                     sequenceColombia, sequenceBoracay, sequenceRioDeJaneiro2, sequenceJamaica, sequenceDominicanRep,
                     sequenceCuiaba, sequenceIndia, sequenceBeloHorizonte, sequenceIndia2, sequenceMexico2,
                     sequenceDNUSBat, sequenceHIVNetherlands, sequenceSepticemia, sequenceXMRV, sequenceHepCEgypt,
                     sequenceHepCFrance, sequenceHIVJapan, sequenceHIVItaly, sequenceDengue, sequenceDengue2,
                     sequenceEncephalitis, sequenceOrf, sequenceOrf2, sequenceHerpes, sequenceSheeppox, sequenceWhiteSpot,
                     sequenceEColi, sequenceNewCastle, sequenceHanta, sequenceIslaVista, sequenceAnemiaVirus,
                     sequenceAnemiaVirus2, sequenceSIV, sequenceSiv2, sequenceHorseHerpes]
    n = len(sequence_list)

    # make a numpy array of size n x n
    upgmaArray = numpy.zeros((n, n))

    # fill array with pairwise distances
    upgmaArray = makeArray(sequence_list, upgmaArray)
    print("UPGMA table for sequences:")
    print(upgmaArray)

    # use the array to get the tree contents
    treeContent = computeTree(upgmaArray)

    # add data to file and draw to tree
    tree = open("tree.dnd", "w+")
    tree.write(treeContent)
    drawTree(treeContent)

    # legend
    print("Phylogenetic Tree Legend")
    print("Original Sequence:")
    print("a-Brazilian Zika Strain Envelope Protein")

    print("\nBasic Data Set")
    print("b-Vietnam Zika Strain Envelope Protein")
    print("c-Puerto Rico Zika Strain Envelope Protein")
    print("d-Guatemala Zika Strain Envelope Protein")
    print("e-Guadeloupe Zika Strain Envelope Protein")
    print("f-USA Zika Strain Envelope Protein")
    print("g-Mexico Zika Strain Envelope Protein")
    print("h-Martinique Zika Strain Envelope Protein")
    print("i-Yucatan Zika Strain Envelope Protein")
    print("j-Haiti Zika Strain Envelope Protein")
    print("k-Fiji Zika Strain Envelope Protein")
    print("l-French Guiana Zika Strain Envelope Protein")
    print("m-Angola Zika Strain Envelope Protein")
    print("n-Cuba Zika Strain Envelope Protein")
    print("o-Honduras Zika Strain Envelope Protein")
    print("p-Rio de Janeiro Zika Strain Envelope Protein")
    print("q-Colombia Zika Strain Envelope Protein")
    print("r-Boracay Zika Strain Envelope Protein")
    print("s-2nd Rio de Janeiro Zika Strain Envelope Protein")
    print("t-Jamaica Zika Strain Envelope Protein")
    print("u-Dominican Republic Zika Strain Envelope Protein")
    print("v-Cuiaba Zika Strain Envelope Protein")
    print("w-India Zika Strain Envelope Protein")
    print("x-Belo Horizonte Zika Strain Envelope Protein")
    print("y-2nd India Zika Strain Envelope Protein")
    print("z-2nd Mexico Zika Strain Envelope Protein")

    print("\nRelated Data Set")
    print("aa-Bat Flavivirus Envelope Protein")
    print("ab-Netherlands HIV Strain Envelope Protein")
    print("ac-Septicemia Envelope Protein")
    print("ad-XMRV Envelope Protein")
    print("ae-Egypt Hepatitis C Strain Envelope Protein")
    print("af-France Hepatitis C Strain Envelope Protein")
    print("ag-Japan HIV Strain Envelope Protein")
    print("ah-Italy HIV Strain Envelope Protein")
    print("ai-Dengue Envelope Protein")
    print("aj-2nd Dengue Envelope Protein")
    print("ak-Encephalitis Envelope Protein")
    print("al-Orf Virus Envelope Protein")
    print("am-2nd Orf Virus Envelope Protein")
    print("an-Herpes Envelope Protein")
    print("ao-Sheeppox Envelope Protein")
    print("ap-Shrimp White Spot Envelope Protein")
    print("aq-E Coli Envelope Protein")
    print("ar-New Castle Virus Envelope Protein")
    print("as-Hanta Virus Envelope Protein")
    print("at-Isla Vista Virus Envelope Protein")
    print("au-Equine Infectious Anemia Envelope Protein")
    print("av-2nd Equine Infectious Anemia Envelope Protein")
    print("aw-SIV Envelope Protein")
    print("ax-2nd SIV Envelope Protein")
    print("ay-Equine Herpes Envelope Protein")


# takes tree file and constructs a tree using biopython
def drawTree(treeFile):
    tree = Phylo.read(StringIO(treeFile), "newick")
    Phylo.draw_ascii(tree)


# makes a upgma array using pairwise distances between sequences
def makeArray(list, array):
    for i, element1 in enumerate(list):
        for j, element2 in enumerate(list):
            if j >= i:
                # since matrix is mirrored, no need to enumerate j past i
                break
            # calculate pairwise distance
            distance = pairwise(element1, element2)
            # add pairwise distances to array
            array[i, j] = distance
            array[j, i] = distance
    return array


# use upgma array to implement upgma algorithm
def computeTree(array):
    # create labels for use in Array and tree
    labels = []
    for i in range(ord("A"), ord("Z") + 1):
        labels.append(chr(i))
    for i in range(ord("A"), ord("Y") + 1):
        labels.append("A" + chr(i))

    # get lower triangular matrix
    arrayData = [
        [],
        [array[1][0]],
        [array[2][0], array[2][1]],
        [array[3][0], array[3][1], array[3][2]],
        [array[4][0], array[4][1], array[4][2], array[4][3]],
        [array[5][0], array[5][1], array[5][2], array[5][3], array[5][4]],
        [array[6][0], array[6][1], array[6][2], array[6][3], array[6][4], array[6][5]],
        [array[7][0], array[7][1], array[7][2], array[7][3], array[7][4], array[7][5], array[7][6]],
        [array[8][0], array[8][1], array[8][2], array[8][3], array[8][4], array[8][5], array[8][6], array[8][7]],
        [array[9][0], array[9][1], array[9][2], array[9][3], array[9][4], array[9][5], array[9][6], array[9][7], array[9][8]],
        [array[10][0], array[10][1], array[10][2], array[10][3], array[10][4], array[10][5], array[10][6], array[10][7], array[10][8], array[10][9]],
        [array[11][0], array[11][1], array[11][2], array[11][3], array[11][4], array[11][5], array[11][6], array[11][7], array[11][8], array[11][9], array[11][10]],
        [array[12][0], array[12][1], array[12][2], array[12][3], array[12][4], array[12][5], array[12][6], array[12][7], array[12][8], array[12][9], array[12][10], array[12][11]],
        [array[13][0], array[13][1], array[13][2], array[13][3], array[13][4], array[13][5], array[13][6], array[13][7], array[13][8], array[13][9], array[13][10], array[13][11], array[13][12]],
        [array[14][0], array[14][1], array[14][2], array[14][3], array[14][4], array[14][5], array[14][6], array[14][7], array[14][8], array[14][9], array[14][10], array[14][11], array[14][12], array[14][13]],
        [array[15][0], array[15][1], array[15][2], array[15][3], array[15][4], array[15][5], array[15][6], array[15][7], array[15][8], array[15][9], array[15][10], array[15][11], array[15][12], array[15][13], array[15][14]],
        [array[16][0], array[16][1], array[16][2], array[16][3], array[16][4], array[16][5], array[16][6], array[16][7], array[16][8], array[16][9], array[16][10], array[16][11], array[16][12], array[16][13], array[16][14], array[16][15]],
        [array[17][0], array[17][1], array[17][2], array[17][3], array[17][4], array[17][5], array[17][6], array[17][7], array[17][8], array[17][9], array[17][10], array[17][11], array[17][12], array[17][13], array[17][14], array[17][15], array[17][16]],
        [array[18][0], array[18][1], array[18][2], array[18][3], array[18][4], array[18][5], array[18][6], array[18][7], array[18][8], array[18][9], array[18][10], array[18][11], array[18][12], array[18][13], array[18][14], array[18][15], array[18][16], array[18][17]],
        [array[19][0], array[19][1], array[19][2], array[19][3], array[19][4], array[19][5], array[19][6], array[19][7], array[19][8], array[19][9], array[19][10], array[19][11], array[19][12], array[19][13], array[19][14], array[19][15], array[19][16], array[19][17], array[19][18]],
        [array[20][0], array[20][1], array[20][2], array[20][3], array[20][4], array[20][5], array[20][6], array[20][7], array[20][8], array[20][9], array[20][10], array[20][11], array[20][12], array[20][13], array[20][14], array[20][15], array[20][16], array[20][17], array[20][18], array[20][19]],
        [array[21][0], array[21][1], array[21][2], array[21][3], array[21][4], array[21][5], array[21][6], array[21][7], array[21][8], array[21][9], array[21][10], array[21][11], array[21][12], array[21][13], array[21][14], array[21][15], array[21][16], array[21][17], array[21][18], array[21][19], array[21][20]],
        [array[22][0], array[22][1], array[22][2], array[22][3], array[22][4], array[22][5], array[22][6], array[22][7], array[22][8], array[22][9], array[22][10], array[22][11], array[22][12], array[22][13], array[22][14], array[22][15], array[22][16], array[22][17], array[22][18], array[22][19], array[22][20], array[22][21]],
        [array[23][0], array[23][1], array[23][2], array[23][3], array[23][4], array[23][5], array[23][6], array[23][7], array[23][8], array[23][9], array[23][10], array[23][11], array[23][12], array[23][13], array[23][14], array[23][15], array[23][16], array[23][17], array[23][18], array[23][19], array[23][20], array[23][21], array[23][22]],
        [array[24][0], array[24][1], array[24][2], array[24][3], array[24][4], array[24][5], array[24][6], array[24][7], array[24][8], array[24][9], array[24][10], array[24][11], array[24][12], array[24][13], array[24][14], array[24][15], array[24][16], array[24][17], array[24][18], array[24][19], array[24][20], array[24][21], array[24][22], array[24][23]],
        [array[25][0], array[25][1], array[25][2], array[25][3], array[25][4], array[25][5], array[25][6], array[25][7], array[25][8], array[25][9], array[25][10], array[25][11], array[25][12], array[25][13], array[25][14], array[25][15], array[25][16], array[25][17], array[25][18], array[25][19], array[25][20], array[25][21], array[25][22], array[25][23], array[25][24]],
        [array[26][0], array[26][1], array[26][2], array[26][3], array[26][4], array[26][5], array[26][6], array[26][7], array[26][8], array[26][9], array[26][10], array[26][11], array[26][12], array[26][13], array[26][14], array[26][15], array[26][16], array[26][17], array[26][18], array[26][19], array[26][20], array[26][21], array[26][22], array[26][23], array[26][24], array[26][25]],
        [array[27][0], array[27][1], array[27][2], array[27][3], array[27][4], array[27][5], array[27][6], array[27][7], array[27][8], array[27][9], array[27][10], array[27][11], array[27][12], array[27][13], array[27][14], array[27][15], array[27][16], array[27][17], array[27][18], array[27][19], array[27][20], array[27][21], array[27][22], array[27][23], array[27][24], array[27][25], array[27][26]],
        [array[28][0], array[28][1], array[28][2], array[28][3], array[28][4], array[28][5], array[28][6], array[28][7], array[28][8], array[28][9], array[28][10], array[28][11], array[28][12], array[28][13], array[28][14], array[28][15], array[28][16], array[28][17], array[28][18], array[28][19], array[28][20], array[28][21], array[28][22], array[28][23], array[28][24], array[28][25], array[28][26], array[28][27]],
        [array[29][0], array[29][1], array[29][2], array[29][3], array[29][4], array[29][5], array[29][6], array[29][7], array[29][8], array[29][9], array[29][10], array[29][11], array[29][12], array[29][13], array[29][14], array[29][15], array[29][16], array[29][17], array[29][18], array[29][19], array[29][20], array[29][21], array[29][22], array[29][23], array[29][24], array[29][25], array[29][26], array[29][27], array[29][28]],
        [array[30][0], array[30][1], array[30][2], array[30][3], array[30][4], array[30][5], array[30][6], array[30][7], array[30][8], array[30][9], array[30][10], array[30][11], array[30][12], array[30][13], array[30][14], array[30][15], array[30][16], array[30][17], array[30][18], array[30][19], array[30][20], array[30][21], array[30][22], array[30][23], array[30][24], array[30][25], array[30][26], array[30][27], array[30][28], array[30][29]],
        [array[31][0], array[31][1], array[31][2], array[31][3], array[31][4], array[31][5], array[31][6], array[31][7], array[31][8], array[31][9], array[31][10], array[31][11], array[31][12], array[31][13], array[31][14], array[31][15], array[31][16], array[31][17], array[31][18], array[31][19], array[31][20], array[31][21], array[31][22], array[31][23], array[31][24], array[31][25], array[31][26], array[31][27], array[31][28], array[31][29], array[31][30]],
        [array[32][0], array[32][1], array[32][2], array[32][3], array[32][4], array[32][5], array[32][6], array[32][7], array[32][8], array[32][9], array[32][10], array[32][11], array[32][12], array[32][13], array[32][14], array[32][15], array[32][16], array[32][17], array[32][18], array[32][19], array[32][20], array[32][21], array[32][22], array[32][23], array[32][24], array[32][25], array[32][26], array[32][27], array[32][28], array[32][29], array[32][30], array[32][31]],
        [array[33][0], array[33][1], array[33][2], array[33][3], array[33][4], array[33][5], array[33][6], array[33][7], array[33][8], array[33][9], array[33][10], array[33][11], array[33][12], array[33][13], array[33][14], array[33][15], array[33][16], array[33][17], array[33][18], array[33][19], array[33][20], array[33][21], array[33][22], array[33][23], array[33][24], array[33][25], array[33][26], array[33][27], array[33][28], array[33][29], array[33][30], array[33][31], array[33][32]],
        [array[34][0], array[34][1], array[34][2], array[34][3], array[34][4], array[34][5], array[34][6], array[34][7], array[34][8], array[34][9], array[34][10], array[34][11], array[34][12], array[34][13], array[34][14], array[34][15], array[34][16], array[34][17], array[34][18], array[34][19], array[34][20], array[34][21], array[34][22], array[34][23], array[34][24], array[34][25], array[34][26], array[34][27], array[34][28], array[34][29], array[34][30], array[34][31], array[34][32], array[34][33]],
        [array[35][0], array[35][1], array[35][2], array[35][3], array[35][4], array[35][5], array[35][6], array[35][7], array[35][8], array[35][9], array[35][10], array[35][11], array[35][12], array[35][13], array[35][14], array[35][15], array[35][16], array[35][17], array[35][18], array[35][19], array[35][20], array[35][21], array[35][22], array[35][23], array[35][24], array[35][25], array[35][26], array[35][27], array[35][28], array[35][29], array[35][30], array[35][31], array[35][32], array[35][33], array[35][34]],
        [array[36][0], array[36][1], array[36][2], array[36][3], array[36][4], array[36][5], array[36][6], array[36][7], array[36][8], array[36][9], array[36][10], array[36][11], array[36][12], array[36][13], array[36][14], array[36][15], array[36][16], array[36][17], array[36][18], array[36][19], array[36][20], array[36][21], array[36][22], array[36][23], array[36][24], array[36][25], array[36][26], array[36][27], array[36][28], array[36][29], array[36][30], array[36][31], array[36][32], array[36][33], array[36][34], array[36][35]],
        [array[37][0], array[37][1], array[37][2], array[37][3], array[37][4], array[37][5], array[37][6], array[37][7], array[37][8], array[37][9], array[37][10], array[37][11], array[37][12], array[37][13], array[37][14], array[37][15], array[37][16], array[37][17], array[37][18], array[37][19], array[37][20], array[37][21], array[37][22], array[37][23], array[37][24], array[37][25], array[37][26], array[37][27], array[37][28], array[37][29], array[37][30], array[37][31], array[37][32], array[37][33], array[37][34], array[37][35], array[37][36]],
        [array[38][0], array[38][1], array[38][2], array[38][3], array[38][4], array[38][5], array[38][6], array[38][7], array[38][8], array[38][9], array[38][10], array[38][11], array[38][12], array[38][13], array[38][14], array[38][15], array[38][16], array[38][17], array[38][18], array[38][19], array[38][20], array[38][21], array[38][22], array[38][23], array[38][24], array[38][25], array[38][26], array[38][27], array[38][28], array[38][29], array[38][30], array[38][31], array[38][32], array[38][33], array[38][34], array[38][35], array[38][36], array[38][37]],
        [array[39][0], array[39][1], array[39][2], array[39][3], array[39][4], array[39][5], array[39][6], array[39][7], array[39][8], array[39][9], array[39][10], array[39][11], array[39][12], array[39][13], array[39][14], array[39][15], array[39][16], array[39][17], array[39][18], array[39][19], array[39][20], array[39][21], array[39][22], array[39][23], array[39][24], array[39][25], array[39][26], array[39][27], array[39][28], array[39][29], array[39][30], array[39][31], array[39][32], array[39][33], array[39][34], array[39][35], array[39][36], array[39][37], array[39][38]],
        [array[40][0], array[40][1], array[40][2], array[40][3], array[40][4], array[40][5], array[40][6], array[40][7], array[40][8], array[40][9], array[40][10], array[40][11], array[40][12], array[40][13], array[40][14], array[40][15], array[40][16], array[40][17], array[40][18], array[40][19], array[40][20], array[40][21], array[40][22], array[40][23], array[40][24], array[40][25], array[40][26], array[40][27], array[40][28], array[40][29], array[40][30], array[40][31], array[40][32], array[40][33], array[40][34], array[40][35], array[40][36], array[40][37], array[40][38], array[40][39]],
        [array[41][0], array[41][1], array[41][2], array[41][3], array[41][4], array[41][5], array[41][6], array[41][7], array[41][8], array[41][9], array[41][10], array[41][11], array[41][12], array[41][13], array[41][14], array[41][15], array[41][16], array[41][17], array[41][18], array[41][19], array[41][20], array[41][21], array[41][22], array[41][23], array[41][24], array[41][25], array[41][26], array[41][27], array[41][28], array[41][29], array[41][30], array[41][31], array[41][32], array[41][33], array[41][34], array[41][35], array[41][36], array[41][37], array[41][38], array[41][39], array[41][40]],
        [array[42][0], array[42][1], array[42][2], array[42][3], array[42][4], array[42][5], array[42][6], array[42][7], array[42][8], array[42][9], array[42][10], array[42][11], array[42][12], array[42][13], array[42][14], array[42][15], array[42][16], array[42][17], array[42][18], array[42][19], array[42][20], array[42][21], array[42][22], array[42][23], array[42][24], array[42][25], array[42][26], array[42][27], array[42][28], array[42][29], array[42][30], array[42][31], array[42][32], array[42][33], array[42][34], array[42][35], array[42][36], array[42][37], array[42][38], array[42][39], array[42][40], array[42][41]],
        [array[43][0], array[43][1], array[43][2], array[43][3], array[43][4], array[43][5], array[43][6], array[43][7], array[43][8], array[43][9], array[43][10], array[43][11], array[43][12], array[43][13], array[43][14], array[43][15], array[43][16], array[43][17], array[43][18], array[43][19], array[43][20], array[43][21], array[43][22], array[43][23], array[43][24], array[43][25], array[43][26], array[43][27], array[43][28], array[43][29], array[43][30], array[43][31], array[43][32], array[43][33], array[43][34], array[43][35], array[43][36], array[43][37], array[43][38], array[43][39], array[43][40], array[43][41], array[43][42]],
        [array[44][0], array[44][1], array[44][2], array[44][3], array[44][4], array[44][5], array[44][6], array[44][7], array[44][8], array[44][9], array[44][10], array[44][11], array[44][12], array[44][13], array[44][14], array[44][15], array[44][16], array[44][17], array[44][18], array[44][19], array[44][20], array[44][21], array[44][22], array[44][23], array[44][24], array[44][25], array[44][26], array[44][27], array[44][28], array[44][29], array[44][30], array[44][31], array[44][32], array[44][33], array[44][34], array[44][35], array[44][36], array[44][37], array[44][38], array[44][39], array[44][40], array[44][41], array[44][42], array[44][43]],
        [array[45][0], array[45][1], array[45][2], array[45][3], array[45][4], array[45][5], array[45][6], array[45][7], array[45][8], array[45][9], array[45][10], array[45][11], array[45][12], array[45][13], array[45][14], array[45][15], array[45][16], array[45][17], array[45][18], array[45][19], array[45][20], array[45][21], array[45][22], array[45][23], array[45][24], array[45][25], array[45][26], array[45][27], array[45][28], array[45][29], array[45][30], array[45][31], array[45][32], array[45][33], array[45][34], array[45][35], array[45][36], array[45][37], array[45][38], array[45][39], array[45][40], array[45][41], array[45][42], array[45][43], array[45][44]],
        [array[46][0], array[46][1], array[46][2], array[46][3], array[46][4], array[46][5], array[46][6], array[46][7], array[46][8], array[46][9], array[46][10], array[46][11], array[46][12], array[46][13], array[46][14], array[46][15], array[46][16], array[46][17], array[46][18], array[46][19], array[46][20], array[46][21], array[46][22], array[46][23], array[46][24], array[46][25], array[46][26], array[46][27], array[46][28], array[46][29], array[46][30], array[46][31], array[46][32], array[46][33], array[46][34], array[46][35], array[46][36], array[46][37], array[46][38], array[46][39], array[46][40], array[46][41], array[46][42], array[46][43], array[46][44], array[46][45]],
        [array[47][0], array[47][1], array[47][2], array[47][3], array[47][4], array[47][5], array[47][6], array[47][7], array[47][8], array[47][9], array[47][10], array[47][11], array[47][12], array[47][13], array[47][14], array[47][15], array[47][16], array[47][17], array[47][18], array[47][19], array[47][20], array[47][21], array[47][22], array[47][23], array[47][24], array[47][25], array[47][26], array[47][27], array[47][28], array[47][29], array[47][30], array[47][31], array[47][32], array[47][33], array[47][34], array[47][35], array[47][36], array[47][37], array[47][38], array[47][39], array[47][40], array[47][41], array[47][42], array[47][43], array[47][44], array[47][45], array[47][46]],
        [array[48][0], array[48][1], array[48][2], array[48][3], array[48][4], array[48][5], array[48][6], array[48][7], array[48][8], array[48][9], array[48][10], array[48][11], array[48][12], array[48][13], array[48][14], array[48][15], array[48][16], array[48][17], array[48][18], array[48][19], array[48][20], array[48][21], array[48][22], array[48][23], array[48][24], array[48][25], array[48][26], array[48][27], array[48][28], array[48][29], array[48][30], array[48][31], array[48][32], array[48][33], array[48][34], array[48][35], array[48][36], array[48][37], array[48][38], array[48][39], array[48][40], array[48][41], array[48][42], array[48][43], array[48][44], array[48][45], array[48][46], array[48][47]],
        [array[49][0], array[49][1], array[49][2], array[49][3], array[49][4], array[49][5], array[49][6], array[49][7], array[49][8], array[49][9], array[49][10], array[49][11], array[49][12], array[49][13], array[49][14], array[49][15], array[49][16], array[49][17], array[49][18], array[49][19], array[49][20], array[49][21], array[49][22], array[49][23], array[49][24], array[49][25], array[49][26], array[49][27], array[49][28], array[49][29], array[49][30], array[49][31], array[49][32], array[49][33], array[49][34], array[49][35], array[49][36], array[49][37], array[49][38], array[49][39], array[49][40], array[49][41], array[49][42], array[49][43], array[49][44], array[49][45], array[49][46], array[49][47], array[49][48]],
        [array[50][0], array[50][1], array[50][2], array[50][3], array[50][4], array[50][5], array[50][6], array[50][7], array[50][8], array[50][9], array[50][10], array[50][11], array[50][12], array[50][13], array[50][14], array[50][15], array[50][16], array[50][17], array[50][18], array[50][19], array[50][20], array[50][21], array[50][22], array[50][23], array[50][24], array[50][25], array[50][26], array[50][27], array[50][28], array[50][29], array[50][30], array[50][31], array[50][32], array[50][33], array[50][34], array[50][35], array[50][36], array[50][37], array[50][38], array[50][39], array[50][40], array[50][41], array[50][42], array[50][43], array[50][44], array[50][45], array[50][46], array[50][47], array[50][48], array[50][49]],
    ]

    # distance array for adding branch lengths
    distances = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0]

    # while labels are available to be joined
    while len(labels) > 1:
        print("\nTable being analyzed:")
        for i in range(len(arrayData)):
            print(arrayData[i])
        # find lowest value
        x, y = lowValue(arrayData)

        # merges array entries for x,y by averaging data
        distancesX = (arrayData[x][y] / 2) - distances[x]
        distancesY = (arrayData[x][y] / 2) - distances[y]
        distances[x] = (arrayData[x][y] / 2)
        distances[y] = (arrayData[x][y] / 2)
        mergeArray(arrayData, x, y, labels, distances, distancesX, distancesY)

    return labels[0]


# returns pairwise distance for compared sequences
def pairwise(seq1, seq2):
    return sum(x != y for x, y in zip(seq1, seq2))


# finds location of lowest value in array
def lowValue(array):
    # initialize lowest value at infinity
    low = float("inf")
    x, y = -1, -1

    # iterate table for lowest value
    for i in range(len(array)):
        for j in range(len(array[i])):
            if array[i][j] < low:
                low = array[i][j]
                x, y = i, j
    return x, y


# merges array entries for x,y by averaging data
def mergeArray(array, x, y, labels, distance, distanceX, distanceY):
    # makes sure array is only worked on left side, swap if not
    if y < x:
        x, y = y, x
        tempY = distanceY
        distanceY = distanceX
        distanceX = tempY

    distanceX = str(distanceX)
    distanceY = str(distanceY)
    # create new label for merged data
    labels[x] = "(" + labels[x] + ":" + distanceX + "," + labels[y] + ":" + distanceY + ")"

    # reconstruct the x row
    row = []
    for i in range(0, x):
        row.append((array[x][i] + array[y][i]) / 2)
    array[x] = row

    # reconstruct x column
    for i in range(x + 1, y):
        array[i][x] = (array[i][x] + array[y][i]) / 2

    for i in range(y + 1, len(array)):
        array[i][x] = (array[i][x] + array[i][y]) / 2
        # since data merged, delete leftover column
        del array[i][y]

    # since data merged, delete leftover row
    del array[y]

    # delete old label
    del labels[y]

    # delete old distance
    del distance[y]


#Runs Phylogenetic Tree Program
if __name__ == '__main__':
    print("Comparing sequences of Zika Strain Envelope Protein...\n")
    sequenceCompare()
