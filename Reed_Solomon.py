import math
import copy
import random

textt="проверка \n  кода \n 0384872735465црыаор..чюфыувпапцушнпашнцупадшнцп"

def nahodim_stepeni(aa):
    inf_polinom1 = []
    for i in range(len(aa)):
        if aa[i] == 0:
            inf_polinom1.append(10000000000)
        else:
            inf_polinom1.append(int(math.log2(aa[i])))
    return inf_polinom1

def slozenie_poly_Gallya(a_jj,stepeni_alfa):
    promez_spepen = []
    promez = []
    for j in range(len(a_jj)):
        for h in range(len(stepeni_alfa[a_jj[j]])):
            promez_spepen.append(stepeni_alfa[a_jj[j]][h])
    for j in promez_spepen:
        if promez_spepen.count(j) % 2 != 0:
            promez.append(j)
    if len(promez) > 0:
        s0 = stepeni_alfa.index(sorted(list(set(promez))))
    else:
        s0=10000000000
    return s0


def pravilo_slozenis(i,q,stepeni_alfa):
    flag = True
    stepen_tek = []
    b = []
    if i < q:
        stepen_tek.append(i)
        return stepen_tek
    elif i == q:
        stepen_tek.append(0)
        stepen_tek.append(1)
        return stepen_tek
    elif i > q:
        for j in range(len(stepeni_alfa[-1])):
            stepen = int(stepeni_alfa[-1][j]) + 1
            if stepen == q:
                stepen_tek.append(0)
                stepen_tek.append(1)
                flag = False
            else:
                stepen_tek.append(stepen)
        for j in stepen_tek:
            if stepen_tek.count(j) % 2 != 0:
                b.append(j)
        return sorted(b)

def Kodirovanie_Reed_Solomon(tex):
    kodovie_slova=[]
    p = 2
    q = 4
    t = 3  # число исправ. ошибок
    m = p ** q
    n = p ** q - 1  #
    k = p ** q - 1 - 2 * t  # Число инф.символов
    r = n - k
    d = r + 1
    tex=' '.join(format(ord(x), 'b') for x in tex).split(" ")
    for i in range(len(tex)):
        if len(tex[i])<3*k:
            tex[i] = tex[i].zfill(3*k)
    for ii in range(len(tex)):
        a=[]
        xax = 0
        for jj in range(0, len(tex[ii]), 3):
            if xax < k+1:
                s = ''.join(tex[ii])
                a.append(2**(int(s[jj: 3+jj], 2)))
            xax += 1
        # a=a+1 #Правило понижения степени
        stepeni_alfa = []
        for i in range(n):
            stepeni_alfa.append(pravilo_slozenis(i,q,stepeni_alfa))
        # Кодирование
        inf_polinom = nahodim_stepeni(a)
        s = []
        for i in range(0, n):
            a_j = []
            for j in range(0, k):
                if inf_polinom[j] != 10000000000:
                    tek_stepen2 = j * i + inf_polinom[j]
                    while tek_stepen2 > n - 1:
                        tek_stepen2 = tek_stepen2 - n
                    a_j.append(tek_stepen2)
            # Сложение в поле Галлуа
            hoho = slozenie_poly_Gallya(a_j, stepeni_alfa)
            if hoho != 10000000000:
                s.append(2 ** (hoho))
            else:
                s.append(0)
        kodovie_slova.append(s)
    return kodovie_slova

def naxodim_ozenky(i,n,v,F_1,stepeni_alfa):
    slozenie_vectora = []
    for j in range(n):
        if v[j] != 10000000000:
            slo_ci = F_1[i][j] + v[j]
            while slo_ci > n - 1:
                slo_ci = slo_ci - n
            slozenie_vectora.append(slo_ci)
    for j in range(len(slozenie_vectora) - 1):
        if slozenie_vectora[j] != 10000000000:
            ny = slozenie_poly_Gallya(slozenie_vectora[j:j + 2], stepeni_alfa)
            slozenie_vectora[j + 1] = ny
        elif slozenie_vectora[j + 1] == 10000000000:
            slozenie_vectora[j + 1] = slozenie_vectora[j]
    return slozenie_vectora


def Dekodirovanie_Reed_Solomon(kodovie_slova1):
    p = 2
    q = 4
    t = 3  # число исправ. ошибок
    m = p ** q
    n = p ** q - 1  #
    k = p ** q - 1 - 2 * t  # Число инф.символов
    r = n - k
    d = r + 1
    te = ""
    # a=a+1 #Правило понижения степени
    stepeni_alfa = []
    for i in range(n):
        stepeni_alfa.append(pravilo_slozenis(i, q, stepeni_alfa))
    F = []
    F_1 = []  # Обратная матрица
    for i in range(n):
        stroka = []
        stroka1 = []
        for j in range(n):
            chislo = i * j
            while chislo > n - 1:
                chislo = chislo - n
            if chislo == 0:
                stroka1.append(chislo)
            else:
                stroka1.append(n - chislo)
            stroka.append(chislo)
        F.append(stroka)
        F_1.append(stroka1)
    for kodii in range(len(kodovie_slova1)):
        s=kodovie_slova1[kodii]
        s_stepeni = nahodim_stepeni(s)
        v = copy.deepcopy(s_stepeni)
        # Нахождение оценки информационного вектора
        a_ozenka = []
        for i in range(n):
            a_ozenka.append(naxodim_ozenky(i,n,v,F_1,stepeni_alfa)[-1])
        # Вычисляем синдром
        c = a_ozenka[k:]
        # Транспонированная проверочная матрица
        H_T = []
        for i in range(len(F_1)):
            H_T.append(F_1[i][k:])
        delta = c[0] - c[1]  # На каком месте ошибка
        if delta == 0:
            ind = random.randrange(k, n, 1)
            oshib = []
            oshib.append(0)
            oshib.append(v[ind])
            if v[ind] == 10000000000:
                v[ind] = 0
            else:
                v[ind] = slozenie_poly_Gallya(oshib, stepeni_alfa)
            # Нахождение оценки информационного вектора
            a_ozenka = []
            for i in range(n):
                a_ozenka.append(naxodim_ozenky(i, n, v, F_1, stepeni_alfa)[-1])
            # Вычисляем синдром
            c = a_ozenka[k:]
            delta = c[0] - c[1]  # На каком месте ошибка
        if delta < 0:
            delta = delta + n
        if delta > 0:
            opredel_vec = H_T[delta]
            counter = 0  # Какая ошибка
            flaag=True
            while opredel_vec != c:
                counter += 1
                for i in range(len(opredel_vec)):
                    ful = opredel_vec[i] + 1
                    while ful > n - 1:
                        ful = ful - n
                    opredel_vec[i] = ful
                if counter>16:
                    print("Вы ввели не подходящую последовательность")
                    flaag=False
                    break
            if flaag==True:
                # Вычиясляем корректор
                corr = []
                for i in range(len(F_1)):
                    ese = F_1[i][delta] + counter
                    while ese > n - 1:
                        ese = ese - n
                    corr.append(ese)
                # Находим информационное слово
                a_inf = []
                for i in range(len(corr)):
                    cet = []
                    if a_ozenka[i] == 10000000000:
                        a_inf.append(corr[i])
                    elif corr[i] == 10000000000:
                        a_inf.append(a_ozenka[i])
                    else:
                        cet.append(a_ozenka[i])
                        cet.append(corr[i])
                        a_inf.append(slozenie_poly_Gallya(cet, stepeni_alfa))
                nazalo_inf = []
                for i in range(k):
                    if a_inf[i] == 10000000000:
                        nazalo_inf.append(0)
                    else:
                        nazalo_inf.append(2 ** a_inf[i])
                obratno = ''
                for i in range(k):
                    kkk = str(bin((int(math.log2(nazalo_inf[i]))))[2:])
                    if len(kkk) < 3:
                        kkk = kkk.zfill(3)
                    obratno += kkk
                te += chr(int(obratno, 2))
    return te

teeer = textt.split('\n')
obshee=[]
for takak in range(len(teeer)):
    print(Kodirovanie_Reed_Solomon(teeer[takak]))
    obshee.append(Kodirovanie_Reed_Solomon(teeer[takak]))

for takak in range(len(teeer)):
    print(Dekodirovanie_Reed_Solomon(obshee[takak]))

lol=[[128, 256, 128, 4096, 2048, 16384, 8192, 256, 128, 0, 16, 8192, 64, 1024, 512], [16384, 16, 128, 16384, 64, 0, 8, 4096, 8, 128, 8, 2048, 32, 2, 16]]
print(Dekodirovanie_Reed_Solomon(lol))