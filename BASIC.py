#!/usr/bin/env python3
"""Compatibility entry point for BASIC 2; prefer the basic-assemble command."""
from basic_assembler.assembly import complement, revcom, splitread, extend, stitch
from basic_assembler.cli import parse_args, config_args, main


if __name__ == '__main__':
    raise SystemExit(main())
