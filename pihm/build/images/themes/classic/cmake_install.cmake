# Install script for directory: /Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/Users/bhattgopal/apps/qgis1.0.2.app/Contents/MacOS")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mAction.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionAddAllToOverview.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionAddIsland.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionAddLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionAddLegend.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionAddOgrLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionAddRasterLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionAddRing.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionAddVertex.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionAddWmsLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionCaptureLine.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionCapturePoint.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionCapturePolygon.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionCheckQgisVersion.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionCollapseTree.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionContextHelp.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionCopySelected.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionCustomProjection.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionDeleteAttribute.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionDeleteSelected.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionDeleteVertex.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionDraw.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionEditCopy.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionEditCut.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionEditPaste.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionExpandTree.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionExportMapServer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionFileExit.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionFileNew.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionFileOpen.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionFilePrint.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionFileSave.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionFileSaveAs.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionFileSmall.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionFolder.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionHelpAbout.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionHelpContents.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionHideAllLayers.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionIdentify.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionInOverview.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionInvertSelection.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionLabel.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionMeasure.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionMeasureArea.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionMoveVertex.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionNewAttribute.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionNewBookmark.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionNewFolder.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionNewVectorLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionOpenTable.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionOptions.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionPan.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionProjectProperties.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionPropertyItem.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionQgisHomePage.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionRemove.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionRemoveAllFromOverview.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionRemoveLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionSaveAsSVG.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionSaveMapAsImage.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionScaleBar.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionSelect.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionSelectedToTop.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionShowAllLayers.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionShowBookmarks.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionShowPluginManager.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionToggleEditing.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionUnselectAttributes.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionZoomFullExtent.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionZoomIn.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionZoomLast.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionZoomOut.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionZoomToLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mActionZoomToSelected.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mIconEditable.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mIconGeometryLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mIconLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mIconLineLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mIconNoPyramid.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mIconPointLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mIconPolygonLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mIconProjectionDisabled.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mIconProjectionEnabled.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mIconProjectionProblem.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mIconProperties.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mIconPyramid.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mIconSymbology.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mIconUnknownLayerType.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mIconWaitingForLayerType.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/classic" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/classic/mMapserverExport.png")