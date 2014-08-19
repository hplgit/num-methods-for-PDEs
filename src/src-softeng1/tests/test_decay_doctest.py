import doctest
import nose.tools as nt
import decay_mod_doctest

def test_decay_mod_doctest():
    failure_count, test_count = doctest.testmod(m=decay_mod_doctest)
    nt.assert_equal(failure_count, 0,
                    msg='%d tests out of %d failed' %
                    (failure_count, test_count))

