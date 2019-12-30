.. include:: links.rst

.. _development:

PypeIt Development Procedures and Guidelines
============================================

We encourage anyone to help us develop the PypeIt code base to better
suit your needs and to improve its algorithms.  If you do so, please
follow this list of procedures and guidelines.  In particular, please
note our :ref:`codeconduct`.

Installation
------------

If you plan to develop/alter the PypeIt code directly, you should
install the software via git.  You can fork the GitHub `repo`_, develop
within your own fork, and then perform a pull request from that fork.
Or you can clone the `repo`_ and work with it directly:

.. code-block:: bash

    git clone https://github.com/pypeit/PypeIt
    cd PypeIt
    python3 setup.py develop

It's important that you use the ``develop`` option when you run the
setup script so that changes made to your local installation are
immediately available when you re-execute the code.

This isn't required, but to simplify some of the commands below, I
assume that you have an environmental variable that points to your
installation of the PypeIt repository.  E.g., add the following to your
`.bash_profile`:

.. code-block:: bash

    export PYPEIT_DIR=$HOME/Work/packages/PypeIt

Branches
--------

PypeIt maintains two persistent branches:

 * ``master``: This is the primary stable version of the code.
   Nominally, this is the version of the code that is installed when
   using `pip`_ and its the version that is used to produce latest
   `documentation`_.  Pull requests to this branch are only done before
   tagging a new release of the code or to perform critical bug
   hotfixes.  The release schedule is decided by the core developers.

 * ``develop``: This is the main development version of the code.  It
   should be stable enough to use, but it may contain experimental,
   unsupported code that is work in progress.

When editing the code, please create a new branch stemming from the
``develop`` branch:

.. code-block:: bash

    cd $PYPEIT_DIR
    git checkout develop
    git pull
    git checkout -b my_new_feature

Development Principles and Communication
----------------------------------------

The main thing to keep in mind when developing for PypeIt is that its
primary use is as an end-to-end reduction pipeline.  This has a few
implications that one should keep in mind:

 * By default, the execution of ``run_pypeit`` should continue either
   until a critical error is raised or the reduction is complete.  **No
   direct interaction with the code should be required at any point**,
   unless explicitly requested by the user.  PypeIt does have some
   interactive components, but most of these are executed via separate
   scripts.

 * Any input needed from the user for your feature should be provided by
   a parameter (preferred) or as a command-line argument.

 * When developing and debugging, you may need to interact with the code
   using `pdb`_ or `IPython.embed`_; however, these instances should be
   removed before performing a pull request.

 * The success or failure of any given procedure must be assessed via
   automatically generated quality-assessment figures (preferred) or via
   scripts that interact with the primary output files.

Feature development in PypeIt is unlikely to be fully independent of
other development activities.  Your feature will likely depend on or
influence the outcome of other modules during the data-reduction
process.  This leads to a few important guidelines:

 * Make sure that your branch is always up-to-date with the ``develop``
   branch.  E.g.:

   .. code-block:: bash

        cd $PYPEIT_DIR
        git checkout develop
        git pull
        git checkout my_new_feature
        git merge --no-ff develop

 * Consider the effects of simultaneous development efforts on your work
   and vice versa.  For example, if you're working on a specific module
   of the code that depends on the result/datamodel of the
   wavelength-calibration module, you should communicate this and find
   out if someone else is developing that module and how/if they're
   changing it.  Depending on the scale of those changes, development
   priorities may need to be worked out to minimize merge conflicts and
   the need to immediately rework/refactor new code.

 * When you're nearly ready to submit a PR, communicate this to the team
   so that any merging order can be worked out to ensure a smooth
   process.  Again, because of the interdependent nature of the code
   base, it's difficult to have more than one significant PR open at any
   given time.

Our primary means of **communication** for development is the `PypeIt
developers Slack <https://pypeit.slack.com>`_.  Contact `X Prochaska`_ for
access.

Testing the Code
----------------

PypeIt has two main methods for testing and verifying the code base,
unit tests and a dedicated development suite.

Development Suite
~~~~~~~~~~~~~~~~~

We have compiled a large suite of data from all of PypeIt's supported
instruments that we use to test that PypeIt is successful for *all*
instruments *anytime* new features are developed.  For access to the
shared Google TeamDrive, please contact `Joe Hennawi`_.

To test PypeIt using the data from the Google TeamDrive:

 * Clone the `PypeIt-development-suite`_ repository:

   .. code-block:: bash

        git clone https://github.com/pypeit/PypeIt-development-suite.git

 * Download/sync the TeamDrive to the repository.  The TeamDrive
   ``CALIBS`` and ``RAW_DATA`` directories should be accessible.  For
   syncing, consider using `rclone`_.

   .. warning::

        The ``RAW_DATA`` directory currently contains about 22 Gb of
        data, and running the develop test below produces about 26 Gb of
        reduced data.

 * Run the test suite on the setups designated for development purposes:

   .. code-block:: bash
        
        cd PypeIt-development-suite
        ./pypeit_test develop

   .. warning::
        
        The current script executes 31 tests.  These tests are mostly
        reductions of example data, but they also include fluxing,
        flexure, and coadding tests.  The execution time is system
        dependent, but you should expect it to take approx. 12 hours.

Unit Tests
~~~~~~~~~~

Unit tests are located in the ``$PYPEIT_DIR/pypeit/tests`` directory.  To run
them, make sure you have `pytest`_ installed and then:

.. code-block:: bash

    cd $PYPEIT_DIR/pypeit/tests
    py.test .

.. warning::

    This generates some files that should be ignored.  **Please do not
    add these test files to the repository.**  We're in the process of
    including some automatic clean-up in the testing functions.

Some of the unit tests require the `Development Suite`_ and/or a set of
"cooked" data products with the expected result produced by testing
PypeIt on specific test data.  The unit tests will skip the appropriate
testing function if the `Development Suite`_ path is not defined or if
there is no ``Cooked/`` directory in the relevant path.

For the unit tests to take advantage of the development suite data,
``PYPEIT_DEV`` must be an environmental variable that points to the root
directory with the development suite data.  For example, include the
following line in your `.bash_profile` file in your home directory:

    .. code-block:: bash

        export PYPEIT_DEV=$HOME/Work/packages/PypeIt-development-suite
        
For unit tests that use the "cooked" data, PypeIt must find a directory
called ``$PYPEIT_DEV/Cooked/``.

Workflow
--------

A typical PypeIt development workflow is as follows:

 * Create a new branch stemming from the ``develop`` branch:

   .. code-block:: bash

        cd $PYPEIT_DIR
        git checkout develop
        git pull
        git checkout -b my_new_feature

 * Develop and debug the feature

 * Run the unit tests, fix any failures, add tests that test your new
   feature(s), and/or modify the tests to accommodate your new feature:

   .. code-block:: bash

        cd $PYPEIT_DIR/pypeit/tests
        py.test .

   The tests should be run so that they have access to the `Development
   Suite`_ (so that it can, e.g., test loading data), but this first
   round of tests can/should be run without the "cooked" output (e.g.,
   either delete or move the ``$PYPEIT_DEV/Cooked/`` directory).

 * Run the `Development Suite`_ and fix any failures:

   .. code-block:: bash
        
        cd $PYPEIT_DEV
        ./pypeit_test develop

 * Build the cooked tar file (e.g., replace x.xx.x with some unique
   version for your branch):
 
   .. code-block:: bash
        
        cd $PYPEIT_DEV
        ./build_cooked x.xx.x

 * Rerun and debug the tests (being sure to edit
   ``$PYPEIT_DIR/pypeit/tests/test_cooked.py`` to accept your cooked
   version):

   .. code-block:: bash

        cd $PYPEIT_DIR/pypeit/tests
        py.test .

 * Edit ``$PYPEIT_DIR/CHANGES.rst`` to reflect your key developments and
   update the API `documentation`_.

   .. code-block:: bash

        cd $PYPEIT_DIR
        ./update_docs

 * Make sure all your edits are committed and pushed to the remote
   repository:

   .. code-block:: bash

        cd $PYPEIT_DIR
        git add -u
        git commit -m 'final prep for PR'
        git push

 * `Submit a Pull Request (PR)
   <https://github.com/pypeit/PypeIt/compare>`_. Unless otherwise
   requested, all PRs should be submitted to the ``develop`` branch.

Pull Request Acceptance Requirements
------------------------------------

Once you've submitted a pull request, we'll review your PR and provide
comments on the code.  The minimum requirements for acceptance of a PR
are as follows:

 * If your PR introduces a new instrument that PypeIt is to support for
   the long term, this instrument *must* be added to the `Development
   Suite`_.  That means raw data should be added to the TeamDrive and a
   relevant test should be added to the ``$PYPEIT_DEV/pypeit_test``
   script (via a PR to the `PypeIt-development-suite`_) such that the
   new instrument is included in list of instruments tested by executing
   ``./pypeit_test develop``).

 * The continuous-integration tests performed by TravisCI on the remote
   repository must pass.

 * You have to post a successful report resulting from your execution of
   both the `Unit Tests`_ and the `Development Suite`_.  You should also
   have uploaded the "cooked" data to the TeamDrive; e.g.:

   .. code-block:: bash

        rclone copyto Cooked_pypeit_dev_vx.xx.x.tar.gz gdv:Cooked_pypeit_dev_vx.xx.x.tar.gz

   .. figure:: figures/tests_success.png

        Example posting of successful tests.

   For hotfixes, these tests can be circumvented at the discretion of
   the core developers in the cases where the hotfix is obviously
   correct.

 * All new methods and classes must be at least minimally documented.
   "Minimally documented" means that each method has a docstring that
   gives at least (1) a one sentence description of the purpose of the
   method, (2) a complete list of the required and optional arguments
   and their meaning, (3) a description of the returned objects, if
   there are any.  Documentation is expected to adhere to `Sphinx`_
   syntax; i.e., the docstrings should be `reStructuredText`_.  We
   accept both `Google-format docstrings`_ and `Numpy-format
   docstrings`_, but the former is preferred.

 * The docstrings for any changes to existing methods that were altered
   must have been modified so that they are up-to-date and accurate.

 * Spurious commented code used for debugging or testing is fine, but
   please let us know if you want it to be kept by adding a relevant
   comment, something like ``# TODO: Keep this around for now``, at the
   beginning of the commented block.  Otherwise, we're likely to remove
   the commented code when we come across it.

 * "Unsupported code," that is code that is experimental and still work
   in progress, should be minimized as much as is reasonable.  The
   relevant code block should be clearly marked as experimental or WIP,
   and it should not be executed by the main PypeIt executable,
   ``run_pypeit``.

 * At least two reviewers must accept the code.

Tagging Protocol
----------------

The core development team will regularly tag "release" versions of the
repository.  Tagging a release version of the code is triggered anytime
a hotfix or the development branch of the code is merged into ``master``.
The tagging process is as follows:

 * The intent to tag a new version is announced via a `New Issue
   <https://github.com/pypeit/PypeIt/issues/new>`_ posted to the repo.
   Core developers then have 24 hours to indicate if there's any reason
   to hold off for a minor but critical development that should first be
   merged into ``develop``.  This "grace" period can be circumvented for
   critical hotfixes to ``master``.

 * A `PR <https://github.com/pypeit/PypeIt/compare>`_ is then issued to
   merge ``develop`` or the hotfix into ``master`` and must meet the
   same `Pull Request Acceptance Requirements`_ when merging new
   branches into ``develop``.  For these merges, particular attention
   should be paid to the accuracy of the `documentation`_ and isolation
   of any code that is "unsupported."

 * Once the PR is accepted *but before being merged into master*, the
   code is tagged as follows (uses `bumpversion`_):

    .. code-block:: bash

        cd $PYPEIT_DIR

        # Edit CHANGES.rst to reflect the ending of a development phase
        vi CHANGES.rst
        git add -u
        git commit -m 'tag CHANGES'

        # Increment the version from 'dev' to a release version and tag
        bumpversion release --verbose --tag --commit

        # Push the changes and the new tag
        git push
        git push --tags

 * The PR is accepted and ``develop`` or the hotfix is merged into
   ``master``.  Hotfix branches are deleted, but the ``develop`` branch
   should not be.

 * The tag is released for `pip`_ installation.

 * A branch is created to advance the version of the code to a
   development version string and to update ``develop`` with the new
   master:

    .. code-block:: bash

        cd $PYPEIT_DIR

        # Branch from master
        git checkout master
        git pull
        git checkout -b devup

        # Increment the version to 'dev'
        bumpversion patch --verbose --commit

        # Edit CHANGES.rst to begin the new dev section
        vi CHANGES.rst
        git add -u
        git commit -m 'CHANGES'

 * A quick PR is issued that pulls ``devup`` into ``develop``.  All of
   the `Pull Request Acceptance Requirements`_ should already be
   satisfied, meaning that the PR should be quickly accepted and merged.

----

This document was developed and mutually agreed upon by: Kyle Westfall, J. Xavier Prochaska, Joseph Hennawi

*Last Modified: 20 Dec 2019*


