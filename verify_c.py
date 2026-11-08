import subprocess
import sys
from abelian_groups_sum import count, powerful_walk

binary = sys.argv[1]
for n in [1, 2, 7, 10, 99, 4096, 10**6 + 3, 10**9, 7 * 10**11, 10**13, 123456789012345, 10**18]:
    out = subprocess.run([binary, str(n)], capture_output=True, text=True, check=True).stdout
    value = int(out.split("value=")[1])
    assert value == count(n), n
    if n <= 10**13:
        assert value == powerful_walk(n), n
print("ok")
