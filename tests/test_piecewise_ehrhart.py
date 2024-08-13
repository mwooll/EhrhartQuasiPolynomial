from sage.all import *
from ehrhart_quasi_polynomial import *

if __name__ == "__main__":
    A = Matrix([[-1, 0], [0, -1], [1, 1], [0, 1]])
    sec_fan = secondary_fan(A)

    max_cones, piecewise_qps = piecewise_ehrhart_quasi_polynomial(A)

    for p in sec_fan.rays():
        print(f"point: {p}")
        vertices = create_polytope_from_matrix(A, p).Vrepresentation()
        ehr_poly = ehrhart_quasi_polynomial(vertices)
        print(f"{ehr_poly=}")
    
        for idx, cone in enumerate(max_cones):
            if p in cone:
                print(f"cone idx: {idx}")
                print(piecewise_qps[idx](x0=p[0], x1=p[1], x2=p[2], x3=p[3]))
            continue
        print()