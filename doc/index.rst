.. currentmodule:: pywatershed

####################################################
pywatershed: a hydrologic model in Python
####################################################

Welcome to the `pywatershed` docs!

pywatershed is path for modernizing legacy USGS hydrology models,
in particular PRMS and GSFlow. It is also a place for experimentation with
hydrologic models in this context.

===========
Status
===========

-----------------------
Current Version: 1.0.0
-----------------------
With pywatershed version 1.0.0, we have faithfully reproduced the NHM
configuration of PRMS. For more information see the
`release notes <https://github.com/EC-USGS/pywatershed/releases/tag/1.0.0`_
and the `extended release notes <https://ec-usgs.github.io/pywatershed/2023/11/14/v1-0-0-overview)`_
for version 1.0.0.

-------------------------------
Upcoming development (in 2024)
-------------------------------
The broad goal is to reproduce GSFLOW couping to MODFLOW 6 via its API. This will include
gridded configurations and cascading flows.
We are also working on reservoir representations in pywatershed.


=================
Getting Started
=================

Please note that you can browse the API reference, developer info, and index
using the table of contents on the left. However, the best way to get started
with pywatershed is to dive in to the example notebooks.

| For introductory example notebooks, look in the `examples/ <https://github.com/EC-USGS/pywatershed/tree/main/examples>`_ directory in the repository. Numbered starting at 00, these are meant to be completed in order. Though no notebook outputs are saved in Github, these notebooks can be easily found and run in WholeTale containers (free but sign-up or log-in required):

.. image:: https://raw.githubusercontent.com/whole-tale/wt-design-docs/master/badges/wholetale-explore.svg
   :target: https://dashboard.wholetale.org

| `WholeTale container for latest release (main branch) <https://dashboard.wholetale.org/run/64ae29e8a887f48b9f173678?tab=metadata>`_
| `WholeTale container for develop branch <https://dashboard.wholetale.org/run/64ae25c3a887f48b9f1735c8?tab=metadata>`_

See `README.md <https://github.com/EC-USGS/pywatershed/tree/develop/README.md>`_ for more details on using WholeTale.


========================
Community Engagement
========================

We value your feedback! Please use `discussions <https://github.com/EC-USGS/pywatershed/discussions>`_
or `issues <https://github.com/EC-USGS/pywatershed/issues>`_ on our Github page. You may also suggest
edits to this documentation or open an issue by clicking on the Github Octocat at the top of the page.
For more in-depth contributions, please start by reading over
the `DEVELOPER.md file <https://github.com/EC-USGS/pywatershed/blob/develop/DEVELOPER.md>`_.

Thank you for your interest.


.. toctree::
   :hidden:
   :caption: About

    About <about.rst>

.. toctree::
   :hidden:
   :caption: API Reference

    API Reference <api.rst>
    Atmosphere <api/atmosphere>
    Hydrology <api/hydrology>
    Parameters <api/parameters>
    Base <api/base>
    Utils <api/utils>

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: For developers

   What’s New <whats-new>
   GitHub repository <https://github.com/EC-USGS/pywatershed>

.. toctree::
   :hidden:
   :caption: Index

    Index <genindex>
..
    Module Index <modindex>
