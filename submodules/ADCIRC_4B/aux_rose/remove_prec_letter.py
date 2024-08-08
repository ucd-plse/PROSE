import argparse
import subprocess

descr = """
remove precision specifying letters from floating point literals in fortran codes
"""

parser = argparse.ArgumentParser(description=descr)
parser.add_argument('-f', metavar='fileslist', type=str, required=True,
                    help='fortran files list')
args = parser.parse_args()

def main(fileslist):
    print("transforming files:", fileslist)
    for i in range(0,10):
        print ('{} of 9'.format(i))
        for j in range(0,10):
            subprocess.run('sed -i "s/{}d{}/{}e{}/g" {}'.format(i,j,i,j, args.f), shell=True)
            subprocess.run('sed -i "s/{}D{}/{}E{}/g" {}'.format(i,j,i,j, args.f), shell=True)
        subprocess.run('sed -i "s/\.d{}/.e{}/g" {}'.format(i,i, args.f), shell=True)
        subprocess.run('sed -i "s/\.D{}/.E{}/g" {}'.format(i,i, args.f), shell=True)



if __name__ == "__main__":
    main(args.f)