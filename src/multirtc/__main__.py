import argparse

from multirtc import multirtc
from multirtc.multimetric import ale, point_target, rle


def main():
    global_parser = argparse.ArgumentParser(
        prog='multirtc',
        description='ISCE3-based multi-sensor RTC and cal/val tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = global_parser.add_subparsers(title='command', help='MultiRTC sub-commands')

    multirtc_parser = multirtc.create_parser(subparsers.add_parser('rtc', help=multirtc.__doc__))
    multirtc_parser.set_defaults(func=multirtc.run)

    ale_parser = ale.create_parser(subparsers.add_parser('ale', help=ale.__doc__))
    ale_parser.set_defaults(func=ale.run)

    rle_parser = rle.create_parser(subparsers.add_parser('rle', help=rle.__doc__))
    rle_parser.set_defaults(func=rle.run)

    pt_parser = point_target.create_parser(subparsers.add_parser('pt', help=point_target.__doc__))
    pt_parser.set_defaults(func=point_target.run)

    args = global_parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
