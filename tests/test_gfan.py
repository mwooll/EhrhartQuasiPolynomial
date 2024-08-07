from pathlib import Path
import shutil

import sage.all
from sage.geometry.cone import Cone
from sage.matrix.constructor import Matrix as create_matrix
from sage.modules.free_module_element import free_module_element 

from ehrhart_quasi_polynomial.gfan import (secondary_fan, compute_secondary_fan,
                                           get_names, gfan_secondaryfan,
                                           retrieve_results)

from unittest import TestCase, main


class TestGfan(TestCase):
    def setUp(self):
        Path("test_folder").mkdir(parents=True, exist_ok=True)
        self.dir = "test_folder"
        self.vertices = [[0, 0], [1, 0], [1, 1], [0, 1]]
        
        ray_0 = free_module_element((-1, -1, 1, 0))
        ray_1 = free_module_element((0, 0, 0, 1))
        ray_2 = free_module_element((1, 1, -1, 0))
        self.expected_fan = {
            'AMBIENT_DIM': 4, 'DIM': 4, 'LINEALITY_DIM': 2,
            'RAYS': [ray_0, ray_1, ray_2],
            'N_RAYS': 3, 
            'LINEALITY_SPACE': create_matrix([[1, 0, 1, 0],
                                              [0, 1, 1, 0]]),
            'ORTH_LINEALITY_SPACE': create_matrix([[1, 1, -1, 0],
                                                   [0, 0, 0, 1]]),
            'F_VECTOR': free_module_element((1, 3, 2)),
            'SIMPLICIAL': 1, 'PURE': 1,
            'CONES': [Cone([[0]]),
                      Cone([ray_0]), Cone([ray_1]), Cone([ray_2]),
                      Cone([ray_0, ray_1]), Cone([ray_1, ray_2])],
            'MAXIMAL_CONES': [Cone([ray_0, ray_1]), Cone([ray_1, ray_2])]
                    }

    def tearDown(self):
        shutil.rmtree("test_folder", ignore_errors=True)

    def test_secondary_fan(self):
        A = create_matrix(self.vertices)
        actual = secondary_fan(A)
        self.assertEqual(self.expected_fan, actual)

    def test_compute_secondary_fan(self):
        actual = compute_secondary_fan(self.vertices, True, self.dir)
        self.assertEqual(self.expected_fan, actual)

    def test_get_names(self):
        expected = (self.dir + "/vertices_0", self.dir + "/output_0")
        self.assertEqual(expected, get_names(self.dir))

    def test_gfan_secondaryfan(self):
        input_name, output_name = get_names(self.dir)
        gfan_secondaryfan(self.vertices, input_name, output_name)
        self.assertTrue(Path(output_name).is_file())

    def test_retrieve_results(self):
        input_name, output_name = get_names(self.dir)
        gfan_secondaryfan(self.vertices, input_name, output_name)

        actual = retrieve_results(output_name)
        self.assertEqual(self.expected_fan, actual)

if __name__ == "__main__":
    main()