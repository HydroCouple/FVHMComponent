# Install script for directory: /Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids

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

FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/action.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/contexthelp.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/delete_selected.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/gis_plain_cursor.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/inline_table.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/layers.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionAddAllToOverview.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionAddLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionAddNonDbLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionAddRasterLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionAddWmsLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionCaptureLine.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionCapturePoint.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionCapturePolygon.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionCustomProjection.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionDraw.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionExportMapServer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionFileExit.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionFileNew.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionFileOpen.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionFilePrint.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionFileSave.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionFileSaveAs.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionHelpAbout.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionHelpContents.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionHideAllLayers.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionIdentify.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionInOverview.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionMeasure.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionNewBookmark.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionNewVectorLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionOpenTable.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionOptions.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionPan.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionProjectProperties.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionQgisHomePage.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionRemoveAllFromOverview.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionRemoveLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionSaveMapAsImage.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionSelect.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionShowAllLayers.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionShowBookmarks.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionShowPluginManager.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionZoomFullExtent.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionZoomIn.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionZoomLast.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionZoomOut.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionZoomToLayer.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mActionZoomToSelected.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mIconProjectionDisabled.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mIconProjectionEnabled.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/mIconProjectionProblem.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/over-.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/over.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/over1.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/qgis32x32.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/remove_all_to_overview.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/show_all_layers1.png")
FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/qgis/themes/nkids" TYPE FILE COMPONENTS "Unspecified" FILES "/Users/bhattgopal/dev/cpp/qgis_1.0.2/images/themes/nkids/show_layerReal.png")