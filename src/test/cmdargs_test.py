#!/usr/bin/env python3

import subprocess
import re

# test to check test_args binary, built with cmake, must be run from
# cmake build directory


def test_arg(args, defined_args):
    regex = re.compile("^(\w*)\s+-\s+value:\s+(.*)\s+set:\s+(0|1)$")

    options = []
    for arg, argname, argval in args:
        if args == 'bool' or args == '-b' or args == 'negbool' or args == '-n':
            options += [args]
        else:
            options += [arg, str(argval)]

    ps = subprocess.run(['./test_args'] + options,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True
                        )

    success = True
    if ps.returncode == 0:
        output = []
        for l in ps.stdout.splitlines():
            # print(l)
            m = regex.match(l)
            if m is not None:
                name = m.group(1)
                value = m.group(2)
                is_set = m.group(3)

                in_defined = [x for x in defined_args if x[0] == name]
                in_args = [x for x in args if x[1] == name]
                if is_set == '0':
                    # check it's not in args
                    if len(in_args) > 0:
                        print("{} unset when in args".format(name))
                        success = False
                    # check its default value
                    if len(in_defined) == 1:
                        if in_defined[0][1] != value:

                            print("{}={} not set to default value {}".format(
                                name, value, in_defined[0][1]))
                            success = False
                    else:
                        print("test error")
                        success = False
                elif is_set == '1':
                    # check it's in args
                    if len(in_args) != 1:
                        print("{} not set when in args".format(name))
                        success = False
                    # check its set value
                    if name == 'double':
                        if float(value) != float(in_args[0][2]):
                            print("{}={} wrong value when set from args {}".format(name,
                                                                                   value,
                                                                                   in_args[0][2]
                                                                                   ))
                            success = False
                    elif name == 'int':
                        if int(value) != int(in_args[0][2]):
                            print("{}={} wrong value when set from args {}".format(name,
                                                                                   value,
                                                                                   in_args[0][2]
                                                                                   ))
                            success = False

                    elif name == 'bool':
                        if value == '1':
                            val = 'true'
                        elif value == '0':
                            val = 'false'
                        if val != in_args[0][2]:
                            print("{}={} wrong value when set from args {}".format(name,
                                                                                   value,
                                                                                   in_args[0][2]
                                                                                   ))
                            success = False
                    elif name == 'negbool':
                        if value == '1':
                            val = 'true'
                        elif value == '0':
                            val = 'false'
                        if val != in_args[0][2]:
                            print("{}={} wrong value when set from args {}".format(name,
                                                                                   value,
                                                                                   in_args[0][2]
                                                                                   ))
                            success = False

                    else:
                        if value != in_args[0][2]:
                            print("{}={} wrong value when set from args {}".format(name,
                                                                                   value,
                                                                                   in_args[0][2]
                                                                                   ))
                            success = False
    else:
        print("wrong return code for options ", args)
        print("stdout: ", ps.stdout)
        print("stderr: ", ps.stderr)
        success = False
        return success

    if not success:
        print(74*"*")
        print("FAIL")
        print("input")
        print(args)
        print("output")
        for l in ps.stdout.splitlines():
            print(l)
        print(74*"*")

    return success


# defined arguments in test,
# long form, default value
defined_args = [('int', '-100'),
                ('bool', '0'),
                ('negbool', '1'),
                ('double', '-1.000000e+05'),
                ('string', 'fart')]


success = True
# 4 types, short form
success = test_arg([('-b', 'bool', 'true')], defined_args) and success
success = test_arg([('-n', 'negbool', 'false')], defined_args) and success
success = test_arg([('-i', 'int', '12')], defined_args) and success
success = test_arg([('-d', 'double', 3.14)], defined_args) and success
success = test_arg([('-s', 'string', 'icle')], defined_args) and success
# 4 types, long form
success = test_arg([('--bool', 'bool', 'true')], defined_args) and success
success = test_arg([('--negbool', 'negbool', 'false')],
                   defined_args) and success
success = test_arg([('--int', 'int', '12')], defined_args) and success
success = test_arg([('--double', 'double', 3.14)], defined_args) and success
success = test_arg([('--string', 'string', 'icle')], defined_args) and success

# 4 combined args
success = test_arg([('--bool', 'bool', 'true'),
                    ('--negbool', 'negbool', 'false'),
                    ('--int', 'int', '-1234567890'),
                    ('--double', 'double', 1e4),
                    ('--string', 'string', '../path/to/things/')], defined_args) and success

# 2 combined args
success = test_arg([
    ('-i', 'int', '+20'),
    ('--double', 'double', 15.321E-3)], defined_args) and success

print()
if success:
    print("Test SUCCESS")
else:
    print("Test FAIL")
