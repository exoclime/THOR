#---- code by Urs Schroffenegger -----------------------------------------------

import re
import pathlib


def get_path_matching_regex(root_dir_, pattern):
    """
    Find all files or directories with name matching the regex in directory

    Parameters
    ----------
        root_dir_: path to directory to search
        pattern: regexp pattern to match, as used by `re` module

    Returns
    -------
        files: list of pathlib Path files path matching the pattern
    """
    matched = get_path_matching_regex_with_groups(root_dir_, pattern)
    files = []
    for m in matched:
        files.append(m['path'])

    return files


def get_path_matching_regex_with_groups(root_dir_, pattern):
    """
    Find all files or directories with name matching the regex in directory,
    returning the groups of match.

    Parameters
    ----------
        root_dir_: path to directory to search
        pattern: regexp pattern to match, as used by `re` module

    Returns
    -------
        files: list of dictionnaries for files matching the pattern
               * 'path': pathlib.Path object for matching groups
               * 'groups': re.Match groups for groups in regexp

    Example
    -------
    in directory with

      README.txt
      file_00.txt
      file_01.txt
      file_42.doc
      file_AA.txt

    calling

      get_path_matching_regex_with_groups(".", "file_(\d+).txt")


    returns:


      [{'path': PosixPath('/.../files_00.txt'), 'groups': ('00',)},
       {'path': PosixPath('/.../files_01.txt'), 'groups': ('01',)}]



    """
    search_regex = re.compile(pattern)
    root_dir = pathlib.Path(root_dir_)

    files = []
    for p in sorted(root_dir.glob("*")):
        m = search_regex.fullmatch(p.name)
        if m is not None:
            files.append({'path': p, 'groups': m.groups()})

    return files
