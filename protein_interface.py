Chain1 = "A"
Chain2 = "B"
dist_threshold = 7
inputfile = "4hhb.pdb"

temp = []       #Stores the whole file as a list in temp
a = []          #Stroes the CA atoms of Chain1
b = []          #Stores the CA atoms of Chain2
c = []          #Stores the Euclid Distance which is < Threshold
d = []          #Stores Chain1 CA atoms whose Euclid Distance < Threshold (includes the repeats)
e = []          #Stores Chain2 CA atoms whose Euclid Distance < Threshold (includes the repeats)

with open(inputfile,'r') as fh1:

    for line in fh1.readlines():
        temp.append(line.split()) #Storing the whole file as a list in temp

    for i in range(len(temp)):
        if temp[i][0] == 'ATOM':            #Getting only ATOM out of the temp list
            if temp[i][2] == 'CA':          #Getting only CA atoms
                if temp[i][4] == Chain1:       #Appending only Chain1 in list a
                    a.append(temp[i])
                if temp[i][4] == Chain2:       #Appending only Chain2 in list b
                    b.append(temp[i])

    #Calculate Euclid distance
    for i in range(len(a)):
        for j in range(len(b)):
            euclid_dist = ((float(b[j][6])-float(a[i][6]))**2 + (float(b[j][7])-float(a[i][7]))**2 + (float(b[j][8])-float(a[i][8]))**2)**0.5
            if euclid_dist < dist_threshold:    #Only appending those values lesser than the given threshold into list c
                c.append(a[i][4]+":"+a[i][3]+"("+a[i][5]+")"+" interacts with "+b[j][4]+":"+b[j][3]+"("+b[j][5]+")")
                d.append(a[i])                  #Appending Chain1 CA atoms whose ED < threshold
                e.append(b[j])                  #Appending Chain2 CA atoms whose ED < threshold

    for i in range(len(c)):
        print(c[i])

d1 = []         #Stores only non-repeats from list d
e1 = []         #Stores only non-repeats from list e

f_helix_e = []  #Stores HELIX of Chain1
f_helix_i = []  #Stores HELIX of Chain2
g_sheet_e = []  #Stores SHEET of Chain1
g_sheet_i = []  #Stores SHEEt of Chain2

h_helix_e = []  #Stores HELIX Chain1 interface atoms
h_helix_i = []  #Stores HELIX Chain2 interface atoms
k_sheet_e = []  #Stores SHEET Chain1 interface atoms
k_sheet_i = []  #Stores SHEET Chain2 interface atoms

#Storing non-repeats in list d1
for i in d:
    if i not in d1:
        d1.append(i)

#Storing non-repeats in list e1
for i in e:
    if i not in e1:
        e1.append(i)

for x in range(len(temp)):
    if temp[x][0] == 'HELIX':
        if temp[x][4] == Chain1:
            f_helix_e.append(temp[x])   #List containing helix of Chain1
        if temp[x][4] == Chain2:
            f_helix_i.append(temp[x])   #List containing helix of Chain2
    if temp[x][0] == 'SHEET':
        if temp[x][5] == Chain1:
            g_sheet_e.append(temp[x])   #List containing sheet of Chain1
        if temp[x][5] == Chain2:
            g_sheet_i.append(temp[x])   #List containing sheet of Chain2

#Finding HELIX Chain1
for i in range(len(d1)):
    for j in range(len(f_helix_e)):
        if int(d1[i][5]) in range(int(f_helix_e[j][5]), int(f_helix_e[j][8])+1):
            h_helix_e.append(f_helix_e[j]) #Contains all HELIX Chain1 atoms

#Finding HELIX Chain2
for i in range(len(e1)):
    for j in range(len(f_helix_i)):
        if int(e1[i][5]) in range(int(f_helix_i[j][5]), int(f_helix_i[j][8])+1):
            h_helix_i.append(f_helix_i[j])  #Contains all HELIX Chain2 atoms

#Finding SHEET Chain1
for i in range(len(d1)):
    for j in range(len(g_sheet_e)):
        if int(d1[i][5]) in range(int(g_sheet_e[j][6]), int(g_sheet_e[j][9])+1):
            k_sheet_e.append(g_sheet_e[j])  #Contains all SHEET Chain1 atoms

#Finding SHEET Chain2
for i in range(len(e1)):
    for j in range(len(g_sheet_i)):
        if int(d1[i][5]) in range(int(g_sheet_i[j][6]), int(g_sheet_i[j][9])+1):
            k_sheet_i.append(g_sheet_i[j])  #Contains all SHEET Chain2 atoms
#print(k_sheet_i, len(k_sheet_i))


#Printing the fractions
print(Chain1+" chain")
print(str(len(h_helix_e))+"/"+str(len(d1))+" of the interface amino acids lying on alpha helices.")
print(str(len(k_sheet_e))+"/"+str(len(d1))+" of the interface amino acids lying on beta sheets.")

print("\n"+Chain2+" chain")
print(str(len(h_helix_i))+"/"+str(len(e1))+" of the interface amino acids lying on alpha helices.")
print(str(len(k_sheet_i))+"/"+str(len(e1))+" of the interface amino acids lying on beta sheets.")


#Printing out E-chain closest interface atom to the right
d1.sort(key = lambda x: int(x[5]))  #Sorting Chain1
print("\n"+Chain1+" chain")
for i in range(len(d1)-1):
    print(d1[i][3]+": "+"closest "+d1[i+1][3]+" at distance "+(str(int(d1[i+1][5])-int(d1[i][5]))))

e1.sort(key = lambda x: int(x[5]))  #Sorting Chain2
print("\n"+Chain2+" chain")
for i in range(len(e1)-1):
    print(e1[i][3]+": "+"closest "+e1[i+1][3]+" at distance "+(str(int(e1[i+1][5])-int(e1[i][5]))))
