# module needed for uploading new Incucyte Data to MySQL (Confluency, Condition ID, Metadata)

import glob
import pandas as pd
import numpy as np

#################################################################################################
# functions for creating and updating Incucyte Confluency Dataset

def getIncucyteDataFilenames(PlateNames, ExperimentNames, ProjectNames, Dates, ListOfIncucyteDataFilenames):
    """Collects all filenames containing raw confluency data by using all collected PlateIDs, ExperimentIDs, ProjectIDs and Dates

    Args: PlateNames(all ID of 96 well plate measured), ExperimentNames(ID used to refer to plates testing same conditions), ProjectNames(ID used to refer to plates measured on same day), Dates(date plate was measured), SavedFilenames(empty)

    Returns: list of filenames of raw confluency data in ListOfIncucyteDataFilenames variable"""
    for ExperimentName, ProjectName, PlateName, Date in zip(ExperimentNames, ProjectNames, PlateNames, Dates):
        Path = '/' ExperimentName + '/' + ProjectName + '/' + PlateName + '/'
        IncucyteDataFilename = Path + 'raw_data/' + str(int(Date)) + '_' + PlateName + '.txt'
        ListOfIncucyteDataFilenames.append(IncucyteDataFilename)
    return ListOfIncucyteDataFilenames
    

def loadIncucyteData(ListOfIncucyteDataFilenames, IncucyteData):
    """Opens all raw confluency data files and loads data into IncucyteData

    Args: SavedFilenames(list of filenames of raw confluency data)

    Returns: IncucyteData variable containing raw confluency data per well and timepoint of each plate"""
    for Filename in ListOfIncucyteDataFilenames:

       
        FileContentsDF = pd.read_csv(Filename, sep='\t',decimal=',', na_values = ['Vessel*', 'Metric*', 'Cell Type*', 'Passage*', 'Notes*','Analysis*'],
                              names = ['Date Time','Elapsed','A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12',
                'B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12',
                'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12',
                'D1','D2','D3','D4','D5','D6','D7','D8','D9','D10','D11','D12',
                'E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','E11','E12',
                'F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11','F12',
                'G1','G2','G3','G4','G5','G6','G7','G8','G9','G10','G11','G12',
                'H1','H2','H3','H4','H5','H6','H7','H8','H9','H10','H11','H12'], usecols = ['Elapsed','A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12',
                'B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12',
                'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12',
                'D1','D2','D3','D4','D5','D6','D7','D8','D9','D10','D11','D12',
                'E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','E11','E12',
                'F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11','F12',
                'G1','G2','G3','G4','G5','G6','G7','G8','G9','G10','G11','G12',
                'H1','H2','H3','H4','H5','H6','H7','H8','H9','H10','H11','H12'])
        FileContentsDF.dropna(inplace = True)  

        FileContentsDF = FileContentsDF.rename(columns=FileContentsDF.iloc[0])
        FileContentsDF = FileContentsDF.iloc[pd.RangeIndex(len(FileContentsDF)).drop(0)]

        WellNames = FileContentsDF.columns.tolist()
        WellNames.pop(0)
        ColumnNames = FileContentsDF.columns.tolist()

        # reset index numbering
        FileContentsDF=FileContentsDF.reset_index()

        # remove unnecessary column
        FileContentsDF= FileContentsDF.drop(['index'], axis =1)

        # creating empty dataframe
        AllIncucyteData = pd.DataFrame(columns= ColumnNames)

        # convert each value from comma to decimal
        for ColumnName in ColumnNames:
            for Index, Value in FileContentsDF.items():
                if Index == ColumnName:
                    AllIncucyteData[ColumnName] = Value.str.replace(',', '.')
        FilenameList = Filename.split('/')
        ProjectID = FilenameList[6]
        ProjectIDChunks=ProjectID.split("_")
        ProjectID="_".join(ProjectIDChunks[:2])
        PlateID = FilenameList[7]


        for Index,Row in AllIncucyteData.iterrows(): # row is a series for a given timepoint containing the confluencies across all wells
            Timepoint = AllIncucyteData.loc[Index,'Elapsed']
            for WellID in WellNames:
                Confluency = Row[WellID]
                if Confluency == ' ':
                    Confluency = '0'
                Output=ProjectID,PlateID,Timepoint,WellID,Confluency
                IncucyteData.append(Output)
    return IncucyteData


############################################################################################################################################################################################
# functions for creating and updating Incucyte Condition ID Table

def checkTableExists(TableName,cursor):
    """Checks whether table already exists in MySQL

    Args: TableName(name of table we're interested in finding)

    Returns: True(if table already exists in MySQL) or False(if table doesn't yet exist in MySQL)"""
    cursor.execute("""
        SELECT COUNT(*)
        FROM information_schema.tables
        WHERE table_name = '{0}'
        """.format(TableName.replace('\'', '\'\'')))
    if cursor.fetchone()[0] == 1:
        # cursor.close()
        return True
    # cursor.close()
    return False


def getIncucytePlatelayouts(PlateNames, ExperimentNames, ProjectNames, Dates, ListofIncucytePlatelayouts):
    """Finds all metadata files by using all collected PlateIDs, ExperimentIDs, ProjectIDs and Dates and loads the name of the files into ListofIncucytePlatelayouts

    Args: PlateNames(all ID of 96 well plate measured), ExperimentNames(ID used to refer to plates testing same conditions), ProjectNames(ID used to refer to plates measured on same day), Dates(date plate was measured), SavedFilenames(empty)

    Returns:list of metadata filenames saved in ListofIncucytePlatelayouts"""
    for ExperimentName, ProjectName, PlateName, Date in zip(ExperimentNames, ProjectNames, PlateNames, Dates):
        Path = '/' + ExperimentName + '/' + ProjectName + '/' + PlateName + '/'
        SavedPlatelayout = Path + 'meta_data/' + PlateName + '_PlateLayout.txt'
        ListofIncucytePlatelayouts.append(SavedPlatelayout)
    return ListofIncucytePlatelayouts

def loadIncucyteMetaData(ListofIncucytePlatelayouts, IncucyteMetaData):
    """Opens all metadata files in ListofIncucytePlatelayouts and loads their contents into IncucyteMetaData

    Args: ListofIncucytePlatelayouts(list of metadata filenames)

    Returns:IncucyteMetaData variable containing metadata each plate"""
    for SinglePlatelayout in ListofIncucytePlatelayouts:
        SinglePlatelayoutDF=pd.read_csv(SinglePlatelayout, sep ='\t', encoding='unicode_escape', decimal=",")
        #remove unnecessary columns
        SinglePlatelayoutDF = SinglePlatelayoutDF.drop(['project_ID', 'plate_ID', 'well_ID'], axis =1)
        #fill input_file with data values from plate layout files
        IncucyteMetaData=pd.concat([IncucyteMetaData,SinglePlatelayoutDF])
    return IncucyteMetaData


def removeNonsenseCellLines(IncucyteMetaData,IncucyteMetaDataUpdated):
    """remove nonsense rows where cell_line = None and returns updated IncucyteMetaData (aka: IncucyteMetaDataUpdated)

    Args: IncucyteMetaData(contains all Incucyte metadata)

    Returns:IncucyteMetaDataUpdated with no nonsense cell lines"""
    for Index, Row in IncucyteMetaData.iterrows():
        if Row['cell_line']!='HAP1':
            IncucyteMetaData = IncucyteMetaData.drop(Index)
            IncucyteMetaDataUpdated=pd.concat([IncucyteMetaDataUpdated,IncucyteMetaData])
    return IncucyteMetaDataUpdated

def assignConditionID(IncucyteMetaDataUpdated, SeriesOfConditions):
    """initialize Condition IDs for cases where no existing Incucyte Condition Table exists

    Args: IncucyteMetaDataUpdated(contains all Incucyte metadata with no nonsense cell lines)

    Returns:SeriesOfConditions(series of newly generated condition IDs)"""
    for Index, Row in IncucyteMetaDataUpdated.iterrows():
        ConditionID = getConditionID(Index)
        SeriesOfConditions.at[Index] = ConditionID
    return SeriesOfConditions

def getConditionID(IndexValue):
    """Creates unique Condition ID string using the IndexValue

    Args: IndexValue(last condition number of last ConditionID listed in MySQL Incucyte Condition Table)

    Returns:ConditionID(string of "condition" + number)"""
    ConditionID = "condition_" + str(IndexValue)
    return ConditionID

def createAdditionalConditionIDs(ExistingAndNewConditionsOverlapDF, SeriesOfConditions, IndexValue, k):
    """Creates ConditionID for each row in ExistingAndNewConditionsOverlapDF and fill SeriesOfConditions with the new IDs based on IndexValue

    Args: ExistingAndNewConditionsOverlapDF(contains conditions already listed in MySQL and newly uploaded conditions),SeriesOfConditions(empty series),IndexValue(last condition number of last ConditionID listed in MySQL Incucyte Condition Table), k(iterator)

    Returns:ExistingAndNewConditionsOverlapDF with newly added ConditionIDs to 'condition_ID' column"""
    for Index, Value in ExistingAndNewConditionsOverlapDF['condition_ID'].items():
        if ExistingAndNewConditionsOverlapDF['condition_ID'].at[Index]==None:
            k +=1
            NewIndexValue = IndexValue + k
            ConditionID = getConditionID(NewIndexValue)
            SeriesOfConditions.at[Index] = ConditionID
            ExistingAndNewConditionsOverlapDF['condition_ID'].at[Index] = SeriesOfConditions.at[Index]
        if Value=='condition_[0-IndexValue]':
            pass
    return ExistingAndNewConditionsOverlapDF

############################################################################################################################################################################################
# functions for creating and updating Incucyte Meta Data Table 

def loadIncucyteMetaData2(ListofIncucytePlatelayouts, IncucyteMetaData):
    """Opens all metadata files in ListofIncucytePlatelayouts and loads their contents into IncucyteMetaData; second version of this function; only difference converts number to int for MySQL

    Args: ListofIncucytePlatelayouts(list of metadata filenames)

    Returns:IncucyteMetaData variable containing metadata each plate"""
    for SinglePlatelayout in ListofIncucytePlatelayouts:
        SinglePlatelayoutDF=pd.read_csv(SinglePlatelayout, sep ='\t', encoding='unicode_escape')
        #fill input_file with data values from plate layout files
        IncucyteMetaData=pd.concat([IncucyteMetaData,SinglePlatelayoutDF])
        IncucyteMetaData['cell_number']=IncucyteMetaData['cell_number'].fillna(0).astype(int)
        # remove rows where well was empty
        IncucyteMetaData.dropna(subset=['cell_line_modifications'], inplace=True)
    return IncucyteMetaData

def createMetaDataTable(IncucyteMetaData,ConditionTableContents):
    """Adds ConditionIDs listed in MySQL Incucyte Condition Table as column to IncucyteMetaData

    Args: IncucyteMetaData(contains all Incucyte metadata), ConditionTableContents(contains entire contents from MySQL Incucyte Condition Table)

    Returns:IncucyteMetaData with repsective ConditionIDs for each well on each plate"""
    MetaDataTable = pd.merge(IncucyteMetaData,ConditionTableContents, on=['cell_line', 'cell_line_modifications','cell_number', 'treatment_0', 'concentration_0', 'treatment_timepoint_0', 
                    'treatment_duration_0','treatment_1', 'concentration_1', 'treatment_timepoint_1', 
                    'treatment_duration_1','treatment_2', 'concentration_2', 'treatment_timepoint_2', 
                    'treatment_duration_2','treatment_3', 'concentration_3', 'treatment_timepoint_3', 
                    'treatment_duration_3'], how="inner")
    # drop unnecessary columns
    MetaDataTable = MetaDataTable.drop(['cell_line','cell_line_modifications','cell_number','treatment_0', 'concentration_0', 'treatment_timepoint_0', 
                        'treatment_duration_0','treatment_1', 'concentration_1', 'treatment_timepoint_1', 
                        'treatment_duration_1','treatment_2','concentration_2', 'treatment_timepoint_2', 
                        'treatment_duration_2','treatment_3', 'concentration_3', 'treatment_timepoint_3', 
                        'treatment_duration_3'], axis =1)
    # remove nans and replace with None
    MetaDataTable = MetaDataTable.replace(np.nan,None)
    # now removing irrelevant rows (aka: conditions not found in my selected plate layout files)
    for Index, Value in MetaDataTable['project_ID'].items():
        if MetaDataTable['project_ID'].at[Index]==None:
            MetaDataTable = MetaDataTable.drop(Index)
        else:
            pass
    for Index, Value in MetaDataTable['plate_ID'].items():
        if MetaDataTable['plate_ID'].at[Index]==None:
            MetaDataTable = MetaDataTable.drop(Index)
        else:
            pass
    for Index, Value in MetaDataTable['well_ID'].items():
        if MetaDataTable['well_ID'].at[Index]==None:
            MetaDataTable = MetaDataTable.drop(Index)
        else:
            pass
    for Index, Value in MetaDataTable['condition_ID'].items():
        if MetaDataTable['condition_ID'].at[Index]==None:
            MetaDataTable = MetaDataTable.drop(Index)
        else:
            pass
    # reset index numbering
    MetaDataTable=MetaDataTable.reset_index()
    # remove unnecessary column
    MetaDataTable = MetaDataTable.drop(['index'], axis =1)
    return MetaDataTable

