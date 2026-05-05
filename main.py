print('Welcome to the tree builder')
print('Please enter a number to select from the options below: ')
print('1 - Create tree with distance matrix from text file')
print('2 - Create tree with distance matrix from user input')
print('3 - Create tree with PHYLIP type multiple sequence alignment file')
dm_option = input('Selection: ')

while dm_option not in ('1', '2', '3'):
    dm_option = input('Please make a selection 1 or 2: ')
    print(dm_option)

if dm_option == '1':
    print('file')
    #distance matrix file feature
    pass

if dm_option == '2':
    print('distance matrix input')
    #Distance matrix input feature
    pass

if dm_option == '3':
    print('PHYLIP style file input')
    #Phylip matrix input feature
    pass

