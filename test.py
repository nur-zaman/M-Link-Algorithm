f = open('data.txt', 'r')
g = open('data2.txt', 'w')
for i in f.readlines():
    value = i.replace('\n', '')
    valueTuple = value.split(' ')
    print(valueTuple)
    # x = f'{int(valueTuple[0])-9} {int(valueTuple[1])-9}\n'
    # print(x)
    g.write(f'{int(valueTuple[0])-9} {int(valueTuple[1])-9}\n')
    # break
