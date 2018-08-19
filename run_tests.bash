set -x

./a.out -b 16 -e 11 > log/b16_e11.txt 2> log/b16_e11_stderr.txt
./a.out -b 16 -e 12 > log/b16_e12.txt 2> log/b16_e12_stderr.txt
./a.out -b 16 -e 13 > log/b16_e13.txt 2> log/b16_e13_stderr.txt
./a.out -b 16 -e 14 > log/b16_e14.txt 2> log/b16_e14_stderr.txt
./a.out -b 16 -e 15 > log/b16_e15.txt 2> log/b16_e15_stderr.txt

./a.out -b 17 -e 12 > log/b17_e12.txt 2> log/b17_e12_stderr.txt
./a.out -b 17 -e 13 > log/b17_e13.txt 2> log/b17_e13_stderr.txt
./a.out -b 17 -e 14 > log/b17_e14.txt 2> log/b17_e14_stderr.txt
./a.out -b 17 -e 15 > log/b17_e15.txt 2> log/b17_e15_stderr.txt
./a.out -b 17 -e 16 > log/b17_e16.txt 2> log/b17_e16_stderr.txt

./a.out -b 18 -e 13 > log/b18_e13.txt 2> log/b18_e13_stderr.txt
./a.out -b 18 -e 14 > log/b18_e14.txt 2> log/b18_e14_stderr.txt
./a.out -b 18 -e 15 > log/b18_e15.txt 2> log/b18_e15_stderr.txt
./a.out -b 18 -e 16 > log/b18_e16.txt 2> log/b18_e16_stderr.txt
./a.out -b 18 -e 17 > log/b18_e17.txt 2> log/b18_e17_stderr.txt

./a.out -b 19 -e 14 > log/b19_e14.txt 2> log/b19_e14_stderr.txt
./a.out -b 19 -e 15 > log/b19_e15.txt 2> log/b19_e15_stderr.txt
./a.out -b 19 -e 16 > log/b19_e16.txt 2> log/b19_e16_stderr.txt
./a.out -b 19 -e 17 > log/b19_e17.txt 2> log/b19_e17_stderr.txt
./a.out -b 19 -e 18 > log/b19_e18.txt 2> log/b19_e18_stderr.txt

./a.out -b 20 -e 15 > log/b20_e15.txt 2> log/b20_e15_stderr.txt
./a.out -b 20 -e 16 > log/b20_e16.txt 2> log/b20_e16_stderr.txt
./a.out -b 20 -e 17 > log/b20_e17.txt 2> log/b20_e17_stderr.txt
./a.out -b 20 -e 18 > log/b20_e18.txt 2> log/b20_e18_stderr.txt
./a.out -b 20 -e 19 > log/b20_e19.txt 2> log/b20_e19_stderr.txt
