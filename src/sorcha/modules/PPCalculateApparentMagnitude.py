from .PPCalculateApparentMagnitudeInFilter import PPCalculateApparentMagnitudeInFilter
from .PPCalculateSimpleCometaryMagnitude import PPCalculateSimpleCometaryMagnitude
from .PPApplyColourOffsets import PPApplyColourOffsets
import logging


def PPCalculateApparentMagnitude(
    observations,
    phasefunction,
    mainfilter,
    othercolours,
    observing_filters,
    object_type,
    lightcurve=False,
    lightcurve_choice="None",
    verbose=False,
):
    """This function applies the correct colour offset to H for the relevant filter, checks to make sure
    the correct columns are included (with additional functionality for colour-specific phase curves),
    then calculates the apparent magnitude.

    Parameters:
    -----------
    observations (Pandas dataframe): dataframe of observations.

    phasefunction (string): desired phase function model. Options are HG, HG12, HG1G2, linear, H.

    mainfilter (string): the main filter in which H is given and all colour offsets are calculated against.

    othercolours (list of strings): list of colour offsets present in input files.

    observing_filters (list of strings): list of observation filters of interest.

    object_type (string): type of object for cometary activity. Either 'comet' or 'none'.

    lightcurve (boolean): whether lightcurves are applied or not

    lc_choice (string): choice of lightcurve model

    verbose (boolean): True/False trigger for verbosity.

    Returns:
    ----------
    observations (Pandas dataframe): dataframe of observations with calculated magnitude column.
    """

    pplogger = logging.getLogger(__name__)
    verboselog = pplogger.info if verbose else lambda *a, **k: None

    if object_type == "comet":
        verboselog("Calculating cometary magnitude using a simple model and applying colour offset...")

        # calculate coma/tail contribution to the apparent magnitude
        observations = PPCalculateSimpleCometaryMagnitude(observations, mainfilter, othercolours)

    # apply correct colour offset to get H magnitude in observation filter
    # if user is only interested in one filter, we have no colour offsets to apply: assume H is in that filter
    if len(observing_filters) > 1:
        verboselog("Selecting and applying correct colour offset...")

        observations = PPApplyColourOffsets(
            observations, phasefunction, othercolours, observing_filters, mainfilter
        )
    else:
        observations.rename(columns={"H_" + mainfilter: "H_filter"}, inplace=True)

    # calculate main body apparent magnitude in observation filter
    verboselog("Calculating apparent magnitude in filter...")
    observations = PPCalculateApparentMagnitudeInFilter(
        observations, phasefunction, lightcurve_choice=lightcurve_choice
    )

    return observations
