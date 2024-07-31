import subprocess
import glob
from pathlib import Path

import sage.all
from sage.geometry.cone import Cone
from sage.matrix.constructor import Matrix as create_matrix
from sage.modules.free_module_element import free_module_element 

def secondary_fan(A, remove=True, directory="gfan_files"):
    return compute_secondary_fan(A.rows(), remove, directory)

def compute_secondary_fan(vectors, remove, directory):
    input_name, output_name = get_names(directory)
    gfan_secondaryfan(vectors, input_name, output_name)
    fan = retrieve_results(output_name)

    if remove:
        command = f"rm {input_name} {output_name}"
        subprocess.Popen(command, shell=True, start_new_session=True).wait()
    return fan 

def get_names(directory):
    Path(directory).mkdir(parents=True, exist_ok=True)
    pre_path = f"{directory}/" if directory else ""

    number = len(glob.glob(f"{pre_path}*"))
    input_name = pre_path + f"vertices_{number}"
    output_name = pre_path + f"output_{number}"
    return input_name, output_name

def gfan_secondaryfan(vertices, input_name, output_name):
    input_vertices = set([tuple(vertex) for vertex in vertices])
    with open(input_name, "w") as file:
        file.write(str(input_vertices))

    command = f"gfan_secondaryfan <{input_name} >{output_name}"
    subprocess.Popen(command, shell=True, start_new_session=True).wait()
    
    with open(output_name, "a") as output:
        output.write("\n")
    return True

def retrieve_results(file_name):
    results = dict()
    name = "undefined"
    value = []
    with open(file_name, "r") as file:
        for line in file.readlines()[4:]:
            stripped = line.strip()

            if stripped == "":
                if len(value) == 1:
                    if len(value[0]) == 1:
                        value = int(value[0])
                results[name] = value
                value = []

            elif stripped[0].isupper(): # name
                name = stripped

            else:
                if "\t" in stripped:
                    index = stripped.index("\t")
                    stripped = stripped[:index]
                value.append(stripped)

    # convert to list of vectors
    for name in ["RAYS", "F_VECTOR"]:
        if name not in results:
            raise ValueError(f"{name} not in results")

        vecs = [free_module_element([int(val) for val in element.split(" ")])
                for element in results[name]]
        results[name] = vecs
        
        if name in ["F_VECTOR"]:
            results[name] = results[name][0]

    # convert to matrices
    for name in ["LINEALITY_SPACE", "ORTH_LINEALITY_SPACE"]:
        if name not in results:
            raise ValueError(f"{name} not is results")
    
        results[name] = create_matrix([[int(val) for val in element.split(" ")]
                                       for element in results[name]])
    
    # convert to cones
    for name in ["CONES", "MAXIMAL_CONES"]:
        if name not in results:
            raise ValueError(f"{name} not in results")
            
        cones = []
        for element in results[name]:
            if element =="{}":
                cones.append(Cone([[0]]))
                continue

            removed = element[1:-1] # removes "{" and "}"
            splitted = removed.split(" ")
            ray = [results["RAYS"][int(val)] for val in splitted]
            cones.append(Cone(ray))
        results[name] = cones
    return results


if __name__ == "__main__":
    vertices = [[1, 0], [0, 1], [0, -1], [-1, -1]]
    result = compute_secondary_fan(vertices, False, "gfan_files")
    for key in result:
        print(f"{key}: {result[key]}")
