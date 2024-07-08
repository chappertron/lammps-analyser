"""
Removes things that are sphinx extensions to the rst standard, so that pandoc can
convert them to markdown properly

Things that are done:
    -  Appends underscores after ` `link text <uri>`_ `
    - These are needed for pandoc to actually recognise these as links,
    according to the rst spec. The format without trailing underscores is a sphinx
    extension
    - Index lines  (`:: index INDEX`) are extracted into a dictionary and removed from
    the file. This dictionary is currently unused but could potentially be turned into
    full proper index pages, or used for lookup of docs files.
    - ` :: parsed-literal` headings are removed from those indented blocks. They are
    interepreted as literals without this in rst.
    -

"""

from os import PathLike
from pathlib import Path

from typing import List
import re


# Just delete these ones...
INDEX_REGEX = re.compile(r".. index::\s+(.+)")
DOCS_REGEX = re.compile(
    r":doc:`(?P<text>.+)\s+<(?P<link_dest>.+)>`"
)  # Replacing with out :doc: and direct link


# TODO: Don't match already proper links (that have trailing _)
LINK_REGEX_MULTI = re.compile(
    r"`(?P<text>[^`]+)\s+<(?P<link_dest>[^`]+)>`"
)  # Multipline version can span multiple lines

# TODO: clear these in the same pass as part of the link regex
DOCS_SIMPLE_REGEX = re.compile(r":doc:")  # Just delete these
# TODO: Convert these to MD footnotes or references.
REF_REGEX = re.compile(
    r":ref:"
)  # Just delete these TODO: Make this work with MD footnotes


def get_file_indices(file_name: PathLike) -> List[str]:
    """
    Extract the  index header from the files.
    """
    with open(file_name, "r", encoding="utf-8") as file:
        file_text = file.read()

    # TODO: Quit early if found the Syntax section
    indices = []
    for match in re.finditer(INDEX_REGEX, file_text):
        indices.append(match.groups()[0])
    return indices


DOCS_PATH = "../lammps_docs/"

# rst index -> file
index_lookup = {}

for file in Path(DOCS_PATH).glob("*.rst"):
    for index in get_file_indices(file):
        # Just the file name
        index_lookup[index] = file.name.removesuffix(".rst")


MODIFIED_PATH = Path("../lammps_docs_cleaned/")


# Go through all the docs files and peform the substitutions
for file in Path(DOCS_PATH).glob("*.rst"):
    print(file)
    with open(file, "r", encoding="utf-8") as stream:
        original_source = stream.read()

    # TODO: Stop applying this regex if gone past the header section?
    modified = INDEX_REGEX.sub("", original_source)

    def replace_link(match: re.Match) -> str:
        # TODO: give these groups names
        text = match.group("text")
        link = match.group("link_dest")
        # file_path = index_lookup[index]

        # TODO: do this with just a sub command and group expansion

        # Make these anonymous links:
        # https://docutils.sourceforge.io/docs/ref/rst/restructuredtext.html#anonymous-hyperlinks
        return rf"`{text} <{link}>`__"

    # Replace any indexes with relative file links
    # modified = DOCS_REGEX.sub(replace_match, modified)

    # Remove any remaining :doc: tags, such as when split accross \n
    modified = DOCS_SIMPLE_REGEX.sub("", modified)

    # Remove any remaining :ref: tags
    modified = REF_REGEX.sub("", modified)

    # Replace any links with proper ones with trailing_
    modified = LINK_REGEX_MULTI.sub(replace_link, modified)

    # Replace with a literal block that can be converted to markdown
    modified = modified.replace(r".. parsed-literal::", r"::")

    with open(MODIFIED_PATH.joinpath(file.name), "w") as stream:
        stream.write(modified)
