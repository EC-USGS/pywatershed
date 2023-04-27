major = 0
minor = 1
micro = 1
__version__ = f"{major}.{minor}.{micro}"

__pakname__ = "pywatershed"

# edit author dictionary as necessary (
# in order of commits after Bakker and Post
author_dict = {
    "James McCreight": "jmccreight@usgs.gov",
    "Christian D. Langevin": "langevin@usgs.gov",
    "Joseph D. Hughes": "jdhughes@usgs.gov",
}
__author__ = ", ".join(author_dict.keys())
__author_email__ = ", ".join(s for _, s in author_dict.items())
