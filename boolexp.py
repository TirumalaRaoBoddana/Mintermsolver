from flask import Flask,render_template,request
from sympy import *
import copy
def countones(binary):
    """It used to find the number of ones in the given binary string.it resturns the count of ones in the binary string"""
    count=0
    for i in range(0,len(binary),1):
        if(binary[i]=="1"):
            count+=1
    return count
def getBinary(number):
    """it we return the binary number of the given number as the string"""
    binary=""
    while(number!=0):
        binary+=str(number%2)
        number=number//2
    return binary[::-1]
def isDifferByOneBinary(a,b):
    """check whether the two binary numbers are differ by 1 binary digit and returns emoty if they are differ by more
    than one 1 difference otherwise return the binary string replaced by the same bit by _ ie 1010,1011 it will return 101_"""
    if len(a) != len(b):
        return ""
    diff=0
    res=""
    for i in range(0,len(a),1):
        if(a[i]!=b[i]):
            res+="_"
            diff+=1
        else:
            res+=a[i]
    if(diff==1):
        return res
    else:
        return ""
def getUnpairedGroups(result,new_group):
    """it will generate the unpaired groups in the each step."""
    unpaired_groups=dict()
    for i in result:
        for j in i[1]:
            pair=list(j.items())[-1][0]
            compared_pairs=[]
            for k in new_group:
                for l in k[1]:
                    compared_pairs.append(list(l.items())[0][0])
            flag=0
            for i in compared_pairs:
                if(set(pair).issubset(set(i))):
                    flag=1
                    break
            if(flag==0):
                unpaired_groups.update({pair:list(j.items())[-1][1]})
    return unpaired_groups
def getCombinedPair(list1,list2):
    """for combining the two groups i.e group 0 and group1, group1 and group2 etc.."""
    pairlist=[]
    for i in list1:
        for j in list2:
            term1 = list(i.items())
            term2 = list(j.items())
            if(isDifferByOneBinary(term1[0][-1],term2[0][-1])!=""):
                pair=term1[0][0]+term2[0][0]
                binaryrepresentation=isDifferByOneBinary(term1[0][-1],term2[0][-1])
                pairlist.append({pair:binaryrepresentation})
    return (len(pairlist),pairlist)
def convertExpression(prime):
    """It will convert the prime implicants to the equivalent boolean expression .ie 10_1 is ab'd """
    result=[]
    for i in prime:
        string=""
        count=0
        for j in i:
            if(j=="_"):
                count+=1
            else:
                if(j=="0"):
                    string=string+chr(97+count)+"'"
                    count+=1
                else:
                    string+=chr(97+count)
                    count+=1
        result.append(string)
    return result
def generatePrimeImplicants(res,unpaired_groups):
    """it takes the result of the quine mcclusky function and returns the tuple which has the prime implicant expression i.e b^c^,bc,ac and corresponding minterms"""
    prime_implicants=[]
    terms=[]
    for i in res:
        for j in i[1]:
            for k in j.items():
                prime_implicants.append(k[1])
                terms.append(k[0])
    for i in unpaired_groups:
        for j in i:
            terms.append(j)
            prime_implicants.append(i[j])
    result_prime=[]
    result_terms=[]
    for i in range(len(prime_implicants)):
        if(prime_implicants[i] not in result_prime):
            result_prime.append(prime_implicants[i])
            result_terms.append(terms[i])
    result_prime=convertExpression(result_prime)
    result_terms=[tuple(sorted(list(i))) for i in result_terms]
    return (result_prime,result_terms)
def countX(list1):
    """It will count the number of X in a column of the prime implicant table here list 1 is the column of the prime implicant chart"""
    index=0
    count=0
    for i in range(len(list1)):
        if(list1[i]=="X"):
            count+=1
            index=i
    if(count==1):
        return index
    else:
        return None
def removeRowsAndColumns(table,index,i,list1,minterms,answer):
    """after finding the column that has single X this function will remove that column from the prime implicant tabel and returns tha table"""
    res=[]
    table1=dict()
    for j in table:
        if(j in list1[-1][index]):
            continue
        else:
            res.append(j)
    for j in res:
        table1.update({j:table[j]})
    for j in table1:
        table1[j].pop(index)
    answer.append(list1[0].pop(index))
    list1[-1].pop(index)
    return (table1,",".join(minterms),list1,answer)
def getExpression(table):
    """it will generate the boolean expression which will covers all the remaining minterms (ex.p=(p1+p2).(p3+p4))"""
    res_exp=""
    for i in table:
        exp=["("]
        for j in range(len(table[i])):
            if(table[i][j]=="X"):
                exp.append("p"+str(j))
        exp="|".join(exp)
        exp+=")"
        exp=list(exp)
        del exp[1]
        exp="".join(exp)
        res_exp+="&"
        res_exp+=exp
    return res_exp.strip("&")
def simplify(var,exp):
    """it is used to simplify the expression obtained by reduced prime implicant table"""
    res=to_dnf(exp)
    res=str(res)
    res=res.split("|")
    result=[]
    for i in res:
        if(len(i.strip(" "))!=len(res[0].strip(" "))):
            break
        else:
            result.append(i)
    return result
def patrikMethod(table,list1):
    #aftre finding the essential prime implicants the implicants and their minterms(list1)
    if(len(table)==0):
        return
    else:
        expression=getExpression(table)
        variables=""
        count=0
        for i in table:
            count=len(table[i])
            break
        for i in range(count):
            variables+="p{0}".format(str(i))
            variables+=","
        variables=variables.strip(",")
        simplified_expression=simplify(variables,expression)
        return simplified_expression
def generateEssentialPrimeImplicants(n,minterms,dontcares,answer):
    """Generate the prime implicants table and identify the essential prime implicants"""
    table = dict()
    res,unpaired_groups,solution=QuineMccluskyAlgo(n,minterms,dontcares)
    implicants=generatePrimeImplicants(res,unpaired_groups)
    for i in minterms.split(","):
        list1 = []
        for j in range(0, len(implicants[1]), 1):
            if int(i) in implicants[1][j]:
                list1.append("X")
            else:
                list1.append("0")
        table.update({int(i): list1})
    prime_implicant_chart=copy.deepcopy(table)
    prime_implicant_terms=copy.deepcopy(list(implicants))
    reduced_chart=[]
    reduced_terms=[]
    while table:
        flag=0
        for i in table:
            if(table[i].count("X")==1):
                flag=1
                index=table[i].index("X")
                table,minterms,list1,answer=removeRowsAndColumns(table,index,i,implicants,minterms,answer)
                reduced_terms.append(copy.deepcopy(list1))
                reduced_chart.append(copy.deepcopy(table))
                break
        if(flag==0):
            break
    expression=""
    result=[]
    if(len(table)>0):
        "here expression is the list of minterms products"
        expression=patrikMethod(table,list1)
        for i in expression:
            minterm_list=[]
            for j in i:
                if(j.isnumeric()):
                    minterm_list.append(implicants[0][int(j)])
            result.append(minterm_list)
        return (result,solution,prime_implicant_chart,prime_implicant_terms,reduced_chart,reduced_terms,expression,answer)
    else:
        return (result,solution,prime_implicant_chart,prime_implicant_terms,reduced_chart,reduced_terms,expression,answer)
def QuineMccluskyAlgo(n,minterms,dontcares):
    """main implementation of quine mcclusky algorithm"""
    terms=[]
    if(len(dontcares)==0):
        try:
            terms=minterms.split(",")
        except:
            return render_template("error.html",error="minterms should not be empty or string")
    else:
        terms=minterms.split(",")+dontcares.split(",")
    groups=dict()
    for i in sorted(terms):
        if(countones(bin(int(i))[2:]) in groups):
            groups[countones(bin(int(i))[2:])].append({(int(i),):getBinary(int(i)).zfill(n)})
        else:
            groups.update({countones(bin(int(i))[2:]):[{(int(i),):getBinary(int(i)).zfill(n)}]})
    res=[i for i in groups]
    res.sort()
    result=dict()
    solution=[]
    for i in res:
        result.update({i:groups[i]})
    solution.append(result)
    result=list(result.items())
    """for storing the unpaired groups"""
    unpaired_groups=[]
    while True:
        new_group=dict()
        for i in range(0,len(result)-1,1):
            list1=result[i][1]
            list2=result[i+1][1]
            if(getCombinedPair(list1,list2)[0]):
                new_group.update({i:getCombinedPair(list1,list2)[-1]})
        if not new_group:
            break
        """reurns the unpaired groups"""
        if(len(getUnpairedGroups(result,list(new_group.items())))!=0):
            unpaired_groups.append(getUnpairedGroups(result,list(new_group.items())))
        solution.append(new_group)
        result=list(new_group.items()).copy()
    return (result,unpaired_groups,solution)
app=Flask(__name__)
@app.route("/")
def homePage():
    return render_template("index.html")
@app.route("/submit",methods=["POST","GET"])
def submitData():
    answer=[]
    if(request.method=="POST"):
        res=request.form
        try:
            var=int(res["var"])
            minterms=res["minterm"].strip(" , ")
            dontcares=res["dontcare"].strip(" , ")
        except ValueError:
            return render_template("error.html",error="number of varibles should be integer")
        if var <= 0:
            return render_template("error.html", error="Number of variables should be greater than 0")
        if not minterms:
            return render_template("error.html", error="Minterms should not be empty")
        try:
            for i in minterms.split(","):
                if(int(i)>=2**var):
                    return render_template("error.html", error="Minterms should be in the range {} to {} only".format(0,2**var-1))
        except:
            return render_template("error.html", error="some error occured")
        if(len(dontcares)!=0):
            for i in dontcares.split(","):
                if(int(i)>=2**var):
                    return render_template("error.html", error="dontcares should be in the range {} to {} only".format(0,2**var-1))
        for i in minterms.split(","):
            if(i in dontcares.split(",")):
                 return render_template("error.html", error="{} is in both minterms and dontcares that is not possible".format(i))
        else:
            non_essential_prime_implicants,solution,prime_implicant_chart,prime_implicant_terms,reduced_chart,reduced_terms,expression,answer=generateEssentialPrimeImplicants(var,minterms,dontcares,answer)
            chart=[]
            for i in range(len(prime_implicant_terms[0])):
                list1=[]
                list1.append(prime_implicant_terms[1][i])
                for j in prime_implicant_chart:
                    list1.append(prime_implicant_chart[j][i])
                list1.append(prime_implicant_terms[0][i])
                chart.append(list1)
            res=""
            for i in set(prime_implicant_terms[0]).difference(set(reduced_terms[0][0])):
                res+=i
            index=prime_implicant_terms[0].index(res)
            terms=prime_implicant_terms[1][index]
            charts=[]
            for i in range(len(reduced_chart)):
                chart1=[]
                for j in range(len(reduced_terms[i][-1])):
                    list1=[]
                    list1.append(reduced_terms[i][-1][j])
                    for k in reduced_chart[i]:
                        list1.append(reduced_chart[i][k][j])
                    list1.append(reduced_terms[i][0][j])
                    chart1.append(list1)
                charts.append(chart1)
            reduced_minterms=[]
            for i in reduced_chart:
                reduced_minterm=[]
                for j in i:
                    reduced_minterm.append(j)
                reduced_minterms.append(reduced_minterm)
            answer1=[]
            if(len(non_essential_prime_implicants)!=0):
                for i in non_essential_prime_implicants:
                    answer1.append(list(set(answer+i)))
            else:
                answer1.append(list(answer))
            simplified_expression=""
            data=[1,2]
            flag=0
            for i in minterms.split(","):
                if(int(i) not in data):
                    flag=1
                    break
            if(flag==0 and var==2):
                simplified_expression+="A⊕B"
            data=[0,3]
            flag=0
            for i in minterms.split(","):
                if(int(i) not in data):
                    flag=1
                    break
            if(flag==0 and var==2):
                simplified_expression+="A⊙B"
            """For three variable special cases"""
            data=[1,2,4,7]
            flag=0
            for i in minterms.split(","):
                if(int(i) not in data):
                    flag=1
                    break
            if(flag==0 and var==3):
                simplified_expression+="A⊕B⊕C"
            data=[0,2,3,5,6]
            flag=0
            for i in minterms.split(","):
                if(int(i) not in data):
                    flag=1
                    break
            if(flag==0 and var==3):
                simplified_expression+="A⊙B⊙C"    
            return render_template("result.html",row=2**(var//2),column=2**(var-var//2)+1,dontcares=dontcares,minterm=[int(i) for i in minterms.split(",")],ans=answer1,sol=solution,prime_implicant_chart=chart,terms=terms,implicant=res,reduced_chart=charts,reduced_minterms=reduced_minterms,expression=expression,simplified_expression=simplified_expression)
    else:
         return render_template("error.html",error="Form should be correctly filled") 
if(__name__=="__main__"):
    app.run(debug=True)