import sys
def main():
    repe=False
    lista=[]
    
    for line in sys.stdin:
        linea = line.split()
        for num in linea:
            if num in lista:
                print("Repetido")
                sys.exit(1)
            else:
                lista.append(num)
    print("All good")
main()