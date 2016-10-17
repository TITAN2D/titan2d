#
# automake template for gui/ directory
#

# Copy files to $(bindir)

if ENABLE_JAVA_GUI
bin_SCRIPTS = titan_gui.sh

jardir=$(prefix)/lib/titan_java_gui
jar_DATA = titan_gui.jar derby.jar jh.jar


JAVAROOT = $(builddir)/
AM_JAVACFLAGS = -cp .:$(srcdir)/derby.jar:$(srcdir)/jh.jar:.
dist_noinst_JAVA = \
    titan/graphics/TextArea.java \
    titan/graphics/DirectorySelector.java \
    titan/graphics/IValueComponent.java \
    titan/graphics/RadioButtonGroup.java \
    titan/graphics/SpringUtilities.java \
    titan/graphics/TextInput.java \
    titan/graphics/LabelDisplay.java \
    titan/graphics/CheckBoxGroup.java \
    titan/gui/JobDetailsDialog.java \
    titan/gui/TabFileIO.java \
    titan/gui/TabGIS.java \
    titan/gui/TabJobMonitor.java \
    titan/gui/DefaultTitanData.java \
    titan/gui/TabFluxSources.java \
    titan/gui/TabPiles.java \
    titan/gui/TabDischargePlanes.java \
    titan/gui/TabMaterialModelMap.java \
    titan/gui/TitanConstants.java \
    titan/gui/TitanTabList.java \
    titan/gui/TabOutputOptions.java \
    titan/gui/Titan.java \
    titan/gui/TabRunParameters.java \
    titan/gui/RunStyleDialog.java \
    titan/gui/TabJobSubmission.java \
    titan/gui/TitanTab.java \
    titan/io/NameValuePairSimple.java \
    titan/io/INameValuePair.java \
    titan/io/NameValuePairComponent.java \
    titan/io/ParameterListInput.java \
    titan/io/ParameterOutput.java \
    titan/io/GenericReader.java \
    titan/io/ParameterInput.java \
    titan/io/NameValuePairGroup.java \
    titan/io/OutputWriter.java \
    titan/io/InputReader.java \
    titan/io/TitanDBAccess.java \
    titan/io/ParameterListOutput.java \
    titan/jobs/RunMethod.java \
    titan/jobs/JobSubmissionParameters.java \
    titan/jobs/JobSubmissionContainer.java \
    titan/jobs/JobStatus.java \
    titan/jobs/JobDetails.java \
    titan/jobs/JobStatusProcessor.java \
    titan/jobs/TitanRunInput.java \
    titan/jobs/condor/RunMethodCondor.java \
    titan/jobs/condor/JobStatusCondor.java \
    titan/jobs/condor/CondorHistoryProcessor.java \
    titan/jobs/condor/CondorJobStatusProcessor.java \
    titan/jobs/condor/CondorQProcessor.java \
    titan/jobs/local/JobStatusLocal.java \
    titan/jobs/local/RunMethodLocal.java \
    titan/jobs/pbs/RunMethodPBS.java \
    titan/jobs/pbs/QstatProcessor.java \
    titan/jobs/pbs/JobStatusPBS.java \
    titan/jobs/submit/RunMethodSubmit.java \
    titan/jobs/submit/JobStatusSubmit.java \
    titan/options/OptionsManager.java \
    titan/options/OptionsDialog.java 

titan_gui.jar:classnoinst.stamp
	$(JAR) -cfe titan_gui.jar titan.gui.Titan titan
endif
