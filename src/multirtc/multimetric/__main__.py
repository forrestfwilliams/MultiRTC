import argparse

from multirtc.multimetric import ale, point_target, rle


def main():
    global_parser = argparse.ArgumentParser(
        prog='multimetric',
        description='SAR cal/val analysis tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = global_parser.add_subparsers(title='analysis', help='SAR metrics')

    ale_parser = ale.create_parser(subparsers.add_parser('ale', help='Absolute Location Error (ALE) analysis'))
    ale_parser.set_defaults(func=ale.run)

    rle_parser = rle.create_parser(subparsers.add_parser('rle', help='Relative Location Error (RLE) analysis'))
    rle_parser.set_defaults(func=rle.run)

    pt_parser = point_target.create_parser(
        subparsers.add_parser('pt', help='Point target (resolution, PSLR, and ISLR) analysis')
    )
    pt_parser.set_defaults(func=point_target.run)

    args = global_parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
