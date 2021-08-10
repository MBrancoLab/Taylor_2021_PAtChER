import csv
import argparse
import sys


def get_args(argv):
    """
    Parse arguments from command line
    """
    args = argparse.ArgumentParser()
    args.add_argument('-u', required = True,
                      metavar = '<bed>',
                      help = 'Patcher output')
    args.add_argument('-r1', required = True,
                      metavar = '<bed>',
                      help = 'Sim3C R1')
    args.add_argument('-r2', required = True,
                      metavar = '<bed>',
                      help = 'Sim3C R2')
    args.add_argument('-wm', action = 'store_true',
                      help = 'Write matched reads to file (default: off)')
    args.add_argument('-wu', action = 'store_true',
                      help = 'Write unmatched reads to file (default: off)')
    return(args.parse_args())


def match_reads(u_file, r1_file, r2_file, write_m, write_u):
    """
    Match IDs between PAtChER reads and Sim3C reads
    Count how many start at the same position
    """

    u = open(u_file, mode='r')
    r1 = open(r1_file, mode='r')
    r2 = open(r2_file, mode='r')

    if write_m:
        wm = open(u_file.replace('.bed', '_matched.bed'), mode='w')
        wm_csv = csv.writer(wm, delimiter='\t')
    if write_u:
        wu = open(u_file.replace('.bed', '_unmatched.bed'), mode='w')
        wu_csv = csv.writer(wu, delimiter='\t')
        wu_csv.writerow(['p.chr','p.start','p.end','id','s.chr','s.start','s.end','dist'])

    u_csv = csv.reader(u, delimiter='\t')
    r1_csv = csv.reader(r1, delimiter='\t')
    r2_csv = csv.reader(r2, delimiter='\t')

    tcount, mcount = 0, 0
    r1_row = r1_csv.__next__()
    r2_row = r2_csv.__next__()

    for row in u_csv:
        rid = row[3].split('_')
        rnum = rid[1].split(':')[0]
        if rnum == '1': #look in R1 file
            while True:
                if r1_row[3] == row[3]:
                    if r1_row[0] == row[0]:
                        dist = abs(int(row[1])-int(r1_row[1]))
                    else:
                        dist = 1e9
                    break
                r1_row = r1_csv.__next__()
        else: #look in R2 file
            while True:
                if r2_row[3] == row[3]:
                    if r2_row[0] == row[0]:
                        dist = abs(int(row[2])-int(r2_row[2]))
                    else:
                        dist = 1e9
                    break
                r2_row = r2_csv.__next__()
            
        tcount += 1
        if dist == 0:
            mcount += 1
            if write_m:
                wm_csv.writerow(row)
        elif write_u:
            if rnum == '1':
                wu_csv.writerow(row[0:4] + r1_row[0:3] + [dist])
            else:
                wu_csv.writerow(row[0:4] + r2_row[0:3] + [dist])

    print('File: ' + u_file)
    print(f'Total PAtChER reads: {tcount}')
    print(f'Total matched reads: {mcount}')

    u.close()
    r1.close()
    r2.close()
    if write_m:
        wm.close()
    if write_u:
        wu.close()


if __name__ == '__main__':
    args = get_args(sys.argv[1:])
    match_reads(u_file=args.u, r1_file=args.r1, r2_file=args.r2, write_m=args.wm, write_u=args.wu)
