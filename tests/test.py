def test_string_for_int(num):
    print(num, type(num))

    if num is None:
        return 1
    try:
        coeff = float(num)
    except (ValueError, TypeError):
        raise ValueError("Coefficient must be numeric")
     
    if not coeff.is_integer():
       raise ValueError("Invalid coefficient, the coefficient cannot be a float") 

    coeff = int(num)

    if coeff > 2:
        raise ValueError("Invalid coefficient, the coordination compound has more than two metal centers.")
    
    return coeff

n=2
m=float(n)
print(m.is_integer())



ligand = "m-ligand"
print(ligand[0:])
