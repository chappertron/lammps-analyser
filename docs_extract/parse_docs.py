from os import PathLike
from pathlib import Path

from typing import List
import re


# TODO: Things to fix
# - {.interpreted-text role="doc"} - These aren't useful?
# - Remove ::: index _ ::: etc.
# - ::: parsed-literal block. -> Turn into nested lists
# - Turn links into references to the markdown files. make a map betwen the index names and file names.

# Just delete these ones...
INDEX_REGEX = re.compile(r".. index::\s+(.+)")
DOCS_REGEX = re.compile(
    r":doc:`(?P<text>.+)\s+<(?P<link_dest>.+)>`"
)  # Replacing with out :doc: and direct link

DOCS_SIMPLE_REGEX = re.compile(r":doc:")  # Just delete these
REF_REGEX = re.compile(
    r":ref:"
)  # Just delete these TODO: Make this work with MD footnotes

# Convert into a markdown style link? or a link to the appropiate file.


# settings = frontend.get_default_settings(Parser)


def get_file_indices(file_name: PathLike) -> List[str]:
    with open(file_name, "r", encoding="utf-8") as file:
        # parser = Parser()
        # document = utils.new_document(file.name, settings)
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

# Because there is a mistake in the lammps docs.
## This index loopup might not actually be needed for linking to files...
index_lookup["fix_atc"] = index_lookup["fix atc"]


MODIFIED_PATH = Path("../lammps_docs_cleaned/")

# print(index_lookup)


# Go through all the docs files and peform the substitutions
for file in Path(DOCS_PATH).glob("*.rst"):
    print(file)
    with open(file, "r", encoding="utf-8") as stream:
        original_source = stream.read()
    # Operate on the string, then write out

    # TODO: Rather than operate line wise, operate on the whole document
    # This will allow for converting the `{text} <{index}>` format to proper links,
    # across multiple lines.
    with open(MODIFIED_PATH.joinpath(file.name), "w") as stream:
        for line in original_source.splitlines():
            # TODO: Stop applying this regex if gone past the header section
            if match := INDEX_REGEX.match(line):
                print(line)
            line = INDEX_REGEX.sub("", line)

            # docs_directives = DOCS_REGEX.findall(line)

            def replace_match(match: re.Match) -> str:
                text = match.group("text")
                index = match.group("link_dest")
                # file_path = index_lookup[index]

                # TODO: do this with just a sub command and group expansion
                return rf"`{text} <{index}>`"

            # Replace any indexes with relative file links
            line = DOCS_REGEX.sub(replace_match, line)

            # Remove any remaining :doc: tags, such as when split accross \n
            line = DOCS_SIMPLE_REGEX.sub("", line)

            # Remove any remaining :ref: tags
            line = REF_REGEX.sub("", line)

            # Replace with a code block that can be converted to markdown
            line = line.replace(r".. parsed-literal::", r"::")

            stream.write(line)
            stream.write("\n")
