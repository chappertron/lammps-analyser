"""
Removes directives from the LAMMPS docs that are Sphinx extensions to the rst standard,
so that pandoc can convert them to markdown properly

Things that are done:
    -  Appends underscores after ` `link text <uri>`_ `
    - These are needed for pandoc to actually recognise these as links,
    according to the rst spec. The format without trailing underscores is a sphinx
    extension
    - Index lines  (`:: index INDEX`) are extracted into a dictionary and removed from
    the file. This dictionary is currently unused but could potentially be turned into
    full proper index pages, or used for lookup of docs files.
    - ` :: parsed-literal` headings are removed from those indented blocks. They are
    interpreted as literals without this in rst.
    -

"""

from os import PathLike
from pathlib import Path

from typing import Dict, List
import re

DOCS_PATH = "../lammps_docs/"
MODIFIED_PATH = Path("../lammps_docs_cleaned/")

# Just delete these ones...
INDEX_REGEX = re.compile(r".. index::\s+(.+)")
DOCS_REGEX = re.compile(
    r":doc:`(?P<text>.+)\s+<(?P<link_dest>.+)>`"
)  # Replacing with out :doc: and direct link


# TODO: Don't match already proper links (that have trailing _)
LINK_REGEX_MULTI = re.compile(
    r"`(?P<text>[^`]+)\s+<(?P<link_dest>[^`]+)>`"
)  # Multipline version can span multiple lines

# TODO: remove these in the one pass, by making part of the link regex
DOCS_SIMPLE_REGEX = re.compile(r":doc:")  # Just delete these
# TODO: Convert these to MD footnotes or references.
REF_REGEX = re.compile(r":ref:")  # Just delete these


def get_file_indices(file_name: PathLike) -> List[str]:
    """
    Extract the  index header from the files.
    """
    with open(file_name, "r", encoding="utf-8") as file:
        file_text = file.read()

    indices = []
    for match in re.finditer(INDEX_REGEX, file_text):
        indices.append(match.groups()[0])
    return indices


def create_index_file_map() -> Dict[str, str]:
    """
    Find index sections in all doc files in each path.
    Create a map from each index to the file that has that index
    """
    # rst index -> file
    index_lookup = {}
    for file in Path(DOCS_PATH).glob("*.rst"):
        for index in get_file_indices(file):
            # Just the file name
            # TODO: What if these aren't a unique map?
            index_lookup[index] = file.name.removesuffix(".rst")
    return index_lookup


def tidy_file(file_contents: str) -> str:
    """
    Removes several Sphinx specific directives from an .rst file

    These are namely:
        -  Appends underscores after ` `link text <uri>`_ `
        These are needed for pandoc to actually recognise these as links,
        according to the [rst spec](https://docutils.sourceforge.io/docs/ref/rst/restructuredtext.html#anonymous-hyperlinks).
        The format without trailing underscores is a
        sphinx extension.
        - Index lines  (`:: index INDEX`) are extracted into a dictionary and removed
        from the file. This dictionary is currently unused but could potentially be
        turned into full proper index pages, or used for lookup of docs files.
        - ` :: parsed-literal` headings are removed from those indented blocks. They
        are interpreted as literals without this in rst.

    ## Implementation:
    The implementation is simply using a few regexes. Ideally a proper .rst parser
    would be used.
    """

    modified = INDEX_REGEX.sub("", file_contents)

    def replace_link(match: re.Match) -> str:
        text = match.group("text")
        link = match.group("link_dest")
        # file_path = index_lookup[index]

        # Make these anonymous links:
        # https://docutils.sourceforge.io/docs/ref/rst/restructuredtext.html#anonymous-hyperlinks
        return rf"`{text} <{link}>`__"

    # Replace any indexes with relative file links
    # modified = DOCS_REGEX.sub(replace_match, modified)

    # Remove any remaining :doc: tags, such as when split across \n
    modified = DOCS_SIMPLE_REGEX.sub("", modified)

    # Remove any remaining :ref: tags
    modified = REF_REGEX.sub("", modified)

    # Replace any links with proper ones with trailing_
    modified = LINK_REGEX_MULTI.sub(replace_link, modified)

    # Replace with a literal block that can be converted to markdown
    modified = modified.replace(r".. parsed-literal::", r"::")

    return modified


def main():
    index_map = create_index_file_map()

    with open("index_map.txt", "w") as stream:
        for k, v in index_map.items():
            stream.write(f"{k},{v}")
            stream.write("\n")

    # Go through all the docs files and perform the substitutions
    for file in Path(DOCS_PATH).glob("*.rst"):
        print(file)
        with open(file, "r", encoding="utf-8") as stream:
            original_source = stream.read()

        modified = tidy_file(original_source)

        with open(MODIFIED_PATH.joinpath(file.name), "w") as stream:
            stream.write(modified)


if __name__ == "__main__":
    main()
