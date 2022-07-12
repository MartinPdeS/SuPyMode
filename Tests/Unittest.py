import unittest

from Validation.Coupler1x1   import CouplerTestCase as Test1
from Validation.Coupler2x2   import CouplerTestCase as Test2
from Validation.Coupler3x3   import CouplerTestCase as Test3
from Validation.Coupler4x4   import CouplerTestCase as Test4


def suite():
    suite = unittest.TestSuite()

    suite.addTest(Test1('test_steps'))

    suite.addTest(Test2('test_steps'))

    suite.addTest(Test3('test_steps'))

    suite.addTest(Test4('test_steps'))

    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner(failfast=True)
    runner.run(suite())








# -
