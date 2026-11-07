cmake_minimum_required(VERSION 3.24)
project(sketchQmlDesktop LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Only the shipped app target (sketch) and its molecule_model library get -O2 below.
# Test executables are left at the default (no optimization flags, since CMAKE_BUILD_TYPE
# is intentionally left unset here) so the edit-build-test loop used constantly during
# development stays as fast as possible -- -O2 on a file the size of IupacNamer.cpp would
# add real seconds to every single-file rebuild for no benefit to a test binary.

find_package(Qt6 6.6 REQUIRED COMPONENTS Core Gui Qml Quick QuickControls2 Concurrent Svg SvgWidgets Test)
qt_standard_project_setup(REQUIRES 6.6)
qt_policy(SET QTP0004 OLD)

# Central naming source list: IupacNamer + its fused-ring geometry helpers
set(IUPAC_NAMING_SOURCES
    src/app/IupacNamer.h src/app/IupacNamer.cpp
    src/app/NamingTrace.h src/app/NamingTrace.cpp
    src/app/FusedRingOrientation.h src/app/FusedRingOrientation.cpp
    src/app/FusedRingDirectionDetector.h src/app/FusedRingDirectionDetector.cpp
    src/app/NamingTypes.h src/app/GraphTypes.h
    src/app/NamingValidation.h src/app/NamingValidation.cpp
    src/app/RetainedRingTypes.h src/app/RetainedRingTypes.cpp
    src/app/RetainedRingMatcher.h src/app/RetainedRingMatcher.cpp
    src/app/PeripheralNumbering.h src/app/PeripheralNumbering.cpp
    src/app/FusionGraphContext.h src/app/FusionGraphContext.cpp
    src/app/NameAssembly.h src/app/NameAssembly.cpp
    src/app/retained/HydrocarbonFamily.h src/app/retained/HydrocarbonFamily.cpp
    src/app/retained/QuinolineFamily.h src/app/retained/QuinolineFamily.cpp
    src/app/retained/IndoleFamily.h src/app/retained/IndoleFamily.cpp
    src/app/retained/TricyclicNFamily.h src/app/retained/TricyclicNFamily.cpp
    src/app/retained/PhenoxazineFamily.h src/app/retained/PhenoxazineFamily.cpp
    src/app/retained/ChromeneXanthineFamily.h src/app/retained/ChromeneXanthineFamily.cpp
    src/app/retained/PhosphorusArsenicFamily.h src/app/retained/PhosphorusArsenicFamily.cpp
)

set_source_files_properties(Theme.qml PROPERTIES QT_QML_SINGLETON_TYPE TRUE)
set_source_files_properties(StyleSheets.qml PROPERTIES QT_QML_SINGLETON_TYPE TRUE)

qt_add_executable(sketch src/app/main.cpp src/app/IndigoService.h src/app/IndigoService.cpp src/app/ImagoService.h src/app/ImagoService.cpp src/v8_process.h src/v8_process.cpp src/app/FileIO.h src/app/PlacementEngines.h src/app/PlacementEngines.cpp src/app/PlacementPreviewManager.h src/app/PlacementPreviewManager.cpp src/app/AppController.h src/app/AppController.cpp src/app/DocumentManager.h src/app/DocumentManager.cpp src/app/molecule/SettingsService.h src/app/molecule/SettingsService.cpp src/app/render/ExportManager.h src/app/render/ExportManager.cpp)
if (WIN32)
    target_sources(sketch PRIVATE src/app/app_icon.rc)
endif()
# FusedRingOrientation/FusedRingDirectionDetector ARE actively used by IupacNamer.cpp's
# N>=3/N>=4 fused-ring numbering logic via computeAllPreferredOrientations,
# computePeripheralNumberingGeneral, detectFusedRingDirections, etc. They are correctly
# included in the sketch target's source list above and are also exercised by their
# own standalone test executables below (confirmed 2026-10-15).

# Console subsystem (the default here) ships with no embedded Windows
# manifest, so Windows falls back to DPI-Unaware and bitmap-stretches the
# whole window to match the display's scale factor -- this is what reads as
# "foggy" text/icons on any scaled display, vs. a native GUI app like Word
# which is manifest-declared Per-Monitor-V2 DPI aware and renders crisp.
# WIN32_EXECUTABLE triggers Qt's CMake tooling to embed that manifest.
# Confirmed via GetProcessDpiAwareness() against the built binary: returned
# Unaware (0) before this change. Trades away the console window (qDebug()/
# stderr no longer show in a terminal) -- run with a debugger attached, or
# temporarily comment this back out, if console output is needed.
set_target_properties(sketch PROPERTIES WIN32_EXECUTABLE TRUE)

qt_add_qml_module(sketch
    URI Sketch.App
    VERSION 1.0
    QML_FILES
        CommandPalette.qml
        MainWindow.qml
        AppMenus.qml
        MainToolbar.qml
        dialogs/FileDialogs.qml
        dialogs/TaskDialogs.qml
        dialogs/MessageDialogs.qml
        dialogs/SettingsDialogs.qml
        ChemCanvas.qml
        MoleculeLayer.qml
        SelectionLayer.qml
        LabelLayer.qml
        ToolPanel.qml
        PropertyPanel.qml
        ToolOverlay.qml
        Theme.qml
        StyleSheets.qml
        PeriodicTablePopup.qml
        BiopolymerDialog.qml
        BiopolymerSequenceView.qml
        AtomPropertiesDialog.qml
        RGroupPanel.qml
        SdfRecordPicker.qml
        SubstructureSearchPopup.qml

        components/FocusRing.qml
        components/IconCell.qml
        components/ThumbnailGridItem.qml
        components/StatusPlaceholder.qml
        components/PopupHeader.qml
        components/AnchoredPicker.qml
        components/SegmentedControl.qml
        components/AppMenuBar.qml
        components/AppMenuBarItem.qml
        components/MenuItemRow.qml
        components/TaskDialog.qml
        components/PreviewCanvas.qml

        js/Selection.js
        js/ToolLabels.js
        js/ThumbnailPainter.js
    RESOURCES
        fonts/materialdesignicons-webfont.ttf
        icons/app_icon.png
        icons/app_icon_16.png
        icons/app_icon_32.png
        icons/app_icon_48.png
        icons/select.svg
        icons/erase.svg
        icons/single_bond.svg
        icons/double_bond.svg
        icons/triple_bond.svg
        icons/up_bond.svg
        icons/down_bond.svg
        icons/updown_bond.svg
        icons/benzene.svg
        icons/clear.svg
        icons/undo.svg
        icons/redo.svg
        icons/copy.svg
        icons/cut.svg
        icons/paste.svg
        icons/save.svg
        icons/open.svg
        icons/fit.svg
        icons/copy_image.svg
        icons/smiles_in.svg
        icons/biopolymer.svg
        icons/select-fragment.svg
        icons/template-0.svg
        icons/template-2.svg
        icons/template-3.svg
        icons/template-4.svg
        icons/template-5.svg
        icons/charge-plus.svg
        icons/charge-minus.svg
        icons/period-table.svg
        icons/generic-groups.svg
        icons/template-lib.svg
        icons/reaction-arrow-open-angle.svg
        icons/reaction-arrow-multitail.svg
        icons/reaction-plus.svg
        icons/reaction-map.svg
        icons/reaction-arrow-filled-triangle.svg
        icons/reaction-arrow-retrosynthetic-arrow.svg
        icons/reaction-arrow-equilibrium-filled-triangle.svg
        icons/reaction-arrow-equilibrium-filled-half-bow.svg
        icons/reaction-arrow-elliptical-arc-arrow-filled-triangle.svg
        icons/chain.svg
        icons/text.svg
        icons/template-6.svg
        icons/template-7.svg
        icons/transform-rotate.svg
        icons/transform-flip-h.svg
        icons/transform-flip-v.svg
        icons/layout.svg
        icons/arom.svg
        icons/dearom.svg
        icons/explicit-hydrogens.svg
        icons/clean.svg
        icons/analyse.svg
        icons/check.svg
        icons/rgroup-label.svg
        icons/zoom-in.svg
        icons/zoom-out.svg
        icons/file-thumbnail.svg
        icons/hand.svg
        icons/select-lasso.svg
        icons/add-image.svg
        icons/search.svg
)

target_link_libraries(sketch PRIVATE Qt6::Core Qt6::Gui Qt6::Qml Qt6::Quick Qt6::QuickControls2 Qt6::Concurrent Qt6::Svg Qt6::SvgWidgets molecule_model iupac_naming
)

# Shipped-binary runtime perf: -O2, scoped to only this target (test executables stay
# unoptimized for fast rebuilds -- see note above CheckIPOSupported). LTO was tried here
# too (INTERPROCEDURAL_OPTIMIZATION) but this MinGW toolchain's ar.exe can't handle the
# LTO plugin objects ("COFF object format mismatch" / "plugin needed to handle lto
# object") -- dropped rather than fight the toolchain for a "cheap win" task.
target_compile_options(sketch PRIVATE -O2)

target_precompile_headers(sketch PRIVATE
    <QtCore/QtCore>
    <QtGui/QtGui>
    <QtQml/QtQml>
    <QtQuick/QtQuick>
    <QtConcurrent/QtConcurrent>
)


# Indigo C++ integration
target_include_directories(sketch PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app/molecule
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app/render
    ${CMAKE_CURRENT_SOURCE_DIR}/build/vcpkg_installed/x64-mingw-dynamic/include
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo-inchi
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo-renderer
    ${CMAKE_CURRENT_SOURCE_DIR}/imago/include
)

if (MSVC)
    # Official EPAM MSVC build (indigo-libs-windows-x86_64, v1.43.0) — the vendored
    # indigo/lib/indigo.dll is MinGW-built (Itanium-mangled exports) and cannot be
    # linked directly by link.exe. indigo.lib/indigo-inchi.lib here are generated import
    # libraries (dumpbin /exports -> .def -> lib.exe /def:), since EPAM's MSVC zip ships
    # only DLLs.
    set(INDIGO_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/indigo-libs-windows-x86_64/lib/windows-x86_64)
    target_link_libraries(sketch PRIVATE ${INDIGO_LIB_DIR}/indigo.lib ${INDIGO_LIB_DIR}/indigo-inchi.lib ${INDIGO_LIB_DIR}/indigo-renderer.lib)
    set(INDIGO_RUNTIME_FILES
        ${INDIGO_LIB_DIR}/indigo.dll
        ${INDIGO_LIB_DIR}/indigo-inchi.dll
        ${INDIGO_LIB_DIR}/indigo-renderer.dll
        ${INDIGO_LIB_DIR}/msvcp140.dll
        ${INDIGO_LIB_DIR}/vcruntime140.dll
        ${INDIGO_LIB_DIR}/vcruntime140_1.dll
    )
else()
    set(INDIGO_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/indigo/lib)
    set(IMAGO_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/imago/lib)
    target_link_libraries(sketch PRIVATE ${INDIGO_LIB_DIR}/indigo.dll ${INDIGO_LIB_DIR}/indigo-inchi.dll ${INDIGO_LIB_DIR}/indigo-renderer.dll ${IMAGO_LIB_DIR}/imago.dll)
    set(INDIGO_RUNTIME_FILES ${INDIGO_LIB_DIR}/indigo.dll ${INDIGO_LIB_DIR}/indigo-inchi.dll ${INDIGO_LIB_DIR}/indigo-renderer.dll ${IMAGO_LIB_DIR}/imago.dll)
endif()

add_custom_command(TARGET sketch POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${INDIGO_RUNTIME_FILES}
        $<TARGET_FILE_DIR:sketch>
)

add_custom_command(TARGET sketch POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/indigo/data/molecules/basic
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/assets/monomer_library.ket
        ${CMAKE_CURRENT_SOURCE_DIR}/indigo/data/molecules/basic/monomer_library.ket
)

qt_add_executable(iupac_namer_test tests/iupac_namer_test.cpp ${IUPAC_NAMING_SOURCES} src/app/BondStereoPerception.h src/app/BondStereoPerception.cpp)
target_include_directories(iupac_namer_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo-inchi
)
if (MSVC)
    target_link_libraries(iupac_namer_test PRIVATE Qt6::Core ${INDIGO_LIB_DIR}/indigo.lib ${INDIGO_LIB_DIR}/indigo-inchi.lib)
else()
    target_link_libraries(iupac_namer_test PRIVATE ${INDIGO_LIB_DIR}/indigo.dll ${INDIGO_LIB_DIR}/indigo-inchi.dll)
endif()

enable_testing()
add_test(NAME iupac_namer_test COMMAND iupac_namer_test)

# Golden-data regression test (Task A/B)
qt_add_executable(iupac_golden_test tests/iupac_golden_test.cpp tools/iupac/AtomOrderRandomizer.h tools/iupac/AtomOrderRandomizer.cpp ${IUPAC_NAMING_SOURCES} src/app/BondStereoPerception.h src/app/BondStereoPerception.cpp)
target_compile_definitions(iupac_golden_test PRIVATE SKETCH_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(iupac_golden_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
    ${CMAKE_CURRENT_SOURCE_DIR}/tools/iupac
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo-inchi
)
if (MSVC)
    target_link_libraries(iupac_golden_test PRIVATE Qt6::Core ${INDIGO_LIB_DIR}/indigo.lib ${INDIGO_LIB_DIR}/indigo-inchi.lib)
else()
    target_link_libraries(iupac_golden_test PRIVATE Qt6::Core ${INDIGO_LIB_DIR}/indigo.dll ${INDIGO_LIB_DIR}/indigo-inchi.dll)
endif()
add_test(NAME iupac_golden_test COMMAND iupac_golden_test)

# Atom order randomization tool (Task A/B)
qt_add_executable(randomize_atom_order tools/iupac/randomize_atom_order.cpp tools/iupac/AtomOrderRandomizer.h tools/iupac/AtomOrderRandomizer.cpp ${IUPAC_NAMING_SOURCES} src/app/BondStereoPerception.h src/app/BondStereoPerception.cpp)
target_include_directories(randomize_atom_order PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
    ${CMAKE_CURRENT_SOURCE_DIR}/tools/iupac
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo-inchi
)
if (MSVC)
    target_link_libraries(randomize_atom_order PRIVATE Qt6::Core ${INDIGO_LIB_DIR}/indigo.lib ${INDIGO_LIB_DIR}/indigo-inchi.lib)
else()
    target_link_libraries(randomize_atom_order PRIVATE Qt6::Core ${INDIGO_LIB_DIR}/indigo.dll ${INDIGO_LIB_DIR}/indigo-inchi.dll)
endif()

qt_add_executable(smiles_driver_v72 tools/iupac/smiles_driver.cpp ${IUPAC_NAMING_SOURCES} src/app/BondStereoPerception.h src/app/BondStereoPerception.cpp)
target_include_directories(smiles_driver_v72 PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo-inchi
)
if (MSVC)
    target_link_libraries(smiles_driver_v72 PRIVATE Qt6::Core ${INDIGO_LIB_DIR}/indigo.lib ${INDIGO_LIB_DIR}/indigo-inchi.lib)
else()
    target_link_libraries(smiles_driver_v72 PRIVATE Qt6::Core ${INDIGO_LIB_DIR}/indigo.dll ${INDIGO_LIB_DIR}/indigo-inchi.dll)
endif()

qt_add_executable(fused_ring_orientation_test tests/fused_ring_orientation_test.cpp src/app/FusedRingOrientation.h src/app/FusedRingOrientation.cpp)
target_include_directories(fused_ring_orientation_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
)
if (MSVC)
    target_link_libraries(fused_ring_orientation_test PRIVATE Qt6::Core)
endif()
add_test(NAME fused_ring_orientation_test COMMAND fused_ring_orientation_test)

qt_add_executable(fused_ring_direction_detector_test tests/fused_ring_direction_detector_test.cpp src/app/FusedRingDirectionDetector.h src/app/FusedRingDirectionDetector.cpp src/app/FusedRingOrientation.h src/app/FusedRingOrientation.cpp)
target_include_directories(fused_ring_direction_detector_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo
)
if (MSVC)
    target_link_libraries(fused_ring_direction_detector_test PRIVATE Qt6::Core ${INDIGO_LIB_DIR}/indigo.lib)
else()
    target_link_libraries(fused_ring_direction_detector_test PRIVATE ${INDIGO_LIB_DIR}/indigo.dll)
endif()
add_test(NAME fused_ring_direction_detector_test COMMAND fused_ring_direction_detector_test)

qt_add_executable(naming_validation_test tests/naming_validation_test.cpp src/app/NamingValidation.h src/app/NamingValidation.cpp src/app/FusedRingOrientation.h src/app/FusedRingOrientation.cpp src/app/NamingTypes.h)
target_include_directories(naming_validation_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
)
if (MSVC)
    target_link_libraries(naming_validation_test PRIVATE Qt6::Core)
else()
    target_link_libraries(naming_validation_test PRIVATE Qt6::Core)
endif()
add_test(NAME naming_validation_test COMMAND naming_validation_test)

qt_add_executable(retained_ring_registry_test tests/retained_ring_registry_test.cpp ${IUPAC_NAMING_SOURCES} src/app/BondStereoPerception.h src/app/BondStereoPerception.cpp)
target_include_directories(retained_ring_registry_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo-inchi
)
if (MSVC)
    target_link_libraries(retained_ring_registry_test PRIVATE Qt6::Core ${INDIGO_LIB_DIR}/indigo.lib ${INDIGO_LIB_DIR}/indigo-inchi.lib)
else()
    target_link_libraries(retained_ring_registry_test PRIVATE Qt6::Core ${INDIGO_LIB_DIR}/indigo.dll ${INDIGO_LIB_DIR}/indigo-inchi.dll)
endif()
add_test(NAME retained_ring_registry_test COMMAND retained_ring_registry_test)

qt_add_executable(retained_ring_matcher_test tests/retained_ring_matcher_test.cpp ${IUPAC_NAMING_SOURCES} src/app/BondStereoPerception.h src/app/BondStereoPerception.cpp)
target_include_directories(retained_ring_matcher_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo-inchi
)
if (MSVC)
    target_link_libraries(retained_ring_matcher_test PRIVATE Qt6::Core ${INDIGO_LIB_DIR}/indigo.lib ${INDIGO_LIB_DIR}/indigo-inchi.lib)
else()
    target_link_libraries(retained_ring_matcher_test PRIVATE Qt6::Core ${INDIGO_LIB_DIR}/indigo.dll ${INDIGO_LIB_DIR}/indigo-inchi.dll)
endif()
add_test(NAME retained_ring_matcher_test COMMAND retained_ring_matcher_test)

add_library(iupac_naming STATIC ${IUPAC_NAMING_SOURCES})
target_include_directories(iupac_naming PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
    ${CMAKE_CURRENT_SOURCE_DIR}/build/vcpkg_installed/x64-mingw-dynamic/include
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo-inchi
)
if (MSVC)
    target_link_libraries(iupac_naming PUBLIC Qt6::Core molecule_model ${INDIGO_LIB_DIR}/indigo.lib ${INDIGO_LIB_DIR}/indigo-inchi.lib)
else()
    target_link_libraries(iupac_naming PUBLIC Qt6::Core molecule_model ${INDIGO_LIB_DIR}/indigo.dll ${INDIGO_LIB_DIR}/indigo-inchi.dll)
endif()
target_compile_options(iupac_naming PRIVATE -O2)

add_library(molecule_model STATIC
    src/app/molecule/ExtensionData.h
    src/app/molecule/MoleculeSnapshot.h
    src/app/molecule/SelectionState.h
    src/app/molecule/EditCommand.h
    src/app/molecule/EditableMolecule.h
    src/app/molecule/EditableMolecule.cpp
    src/app/molecule/DocumentState.h
    src/app/molecule/DocumentState.cpp
    src/app/molecule/TemplateLibrary.h
    src/app/molecule/TemplateLibrary.cpp
    src/app/molecule/ElementData.h
    src/app/molecule/RenderPrimitives.h
    src/app/molecule/RenderPrimitives.cpp
    src/app/molecule/SdfBatch.h
    src/app/molecule/SdfBatch.cpp
    src/app/molecule/ClipboardPreview.h
    src/app/molecule/ClipboardPreview.cpp
    src/app/molecule/BiopolymerSequenceView.h
    src/app/molecule/BiopolymerSequenceView.cpp
    src/app/molecule/DrawingStyle.h
    src/app/molecule/DrawingStyle.cpp
    src/app/molecule/PageSettings.h
    src/app/molecule/PageSettings.cpp
    src/app/molecule/ProjectFile.h
    src/app/molecule/ProjectFile.cpp
    src/app/render/RenderScene.h
    src/app/render/RenderScene.cpp
    src/app/render/Exporter.h
    src/app/render/Exporter.cpp
    src/app/molecule/RenderPrimitivesToVariant.h
    src/app/molecule/RenderPrimitivesToVariant.cpp
    src/app/BondAngleSuggester.h
    src/app/BondAngleSuggester.cpp
    src/app/BondStereoPerception.h
    src/app/BondStereoPerception.cpp
    src/app/PlacementUtils.h
    src/app/PlacementUtils.cpp
    src/app/PlacementCleanup.h
    src/app/PlacementCleanup.cpp
)
target_include_directories(molecule_model PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo
    ${CMAKE_CURRENT_SOURCE_DIR}/indigo/api/c/indigo-renderer
)
if (MSVC)
    target_link_libraries(molecule_model PUBLIC Qt6::Core Qt6::Gui Qt6::Svg Qt6::Concurrent ${INDIGO_LIB_DIR}/indigo.lib ${INDIGO_LIB_DIR}/indigo-renderer.lib)
else()
    target_link_libraries(molecule_model PUBLIC Qt6::Core Qt6::Gui Qt6::Svg Qt6::Concurrent ${INDIGO_LIB_DIR}/indigo.dll ${INDIGO_LIB_DIR}/indigo-renderer.dll)
endif()

# Fewer, larger TUs to compile (no optimization-level change, so this is a pure build-time
# win with no effect on the test executables that link this library unoptimized).
set_target_properties(molecule_model PROPERTIES UNITY_BUILD ON)

# This library is linked into sketch (the shipped app) as well as several test
# executables below; give it the same -O2 as sketch (LTO dropped, see note on the
# sketch target above -- this MinGW toolchain's ar.exe can't handle it).
target_compile_options(molecule_model PRIVATE -O2)

target_precompile_headers(molecule_model PRIVATE <QtCore/QtCore>)

qt_add_executable(editable_molecule_test tests/editable_molecule_test.cpp)
target_link_libraries(editable_molecule_test PRIVATE molecule_model)
target_compile_definitions(editable_molecule_test PRIVATE SKETCH_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME editable_molecule_test COMMAND editable_molecule_test)

qt_add_executable(placement_engines_test tests/placement_engines_test.cpp src/app/PlacementEngines.h src/app/PlacementEngines.cpp)
target_include_directories(placement_engines_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
)
target_link_libraries(placement_engines_test PRIVATE molecule_model)
add_test(NAME placement_engines_test COMMAND placement_engines_test)

# ProjectFile test
qt_add_executable(project_file_test tests/project_file_test.cpp)
target_include_directories(project_file_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
)
target_link_libraries(project_file_test PRIVATE Qt6::Core Qt6::Gui molecule_model)
add_test(NAME project_file_test COMMAND project_file_test)

# SettingsService test
qt_add_executable(settings_service_test tests/settings_service_test.cpp src/app/molecule/SettingsService.h src/app/molecule/SettingsService.cpp)
target_include_directories(settings_service_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
)
target_link_libraries(settings_service_test PRIVATE Qt6::Core Qt6::Gui Qt6::Test molecule_model)
add_test(NAME settings_service_test COMMAND settings_service_test)

# RenderScene test
qt_add_executable(render_scene_test tests/render_scene_test.cpp)
target_include_directories(render_scene_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
)
target_link_libraries(render_scene_test PRIVATE Qt6::Core Qt6::Gui molecule_model)
add_test(NAME render_scene_test COMMAND render_scene_test)

# ExportManager test
qt_add_executable(export_manager_test tests/export_manager_test.cpp src/app/render/ExportManager.h src/app/render/ExportManager.cpp)
target_include_directories(export_manager_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
)
target_link_libraries(export_manager_test PRIVATE Qt6::Core Qt6::Gui Qt6::Svg Qt6::Test molecule_model)
add_test(NAME export_manager_test COMMAND export_manager_test)

qt_add_executable(bond_angle_suggester_test tests/bond_angle_suggester_test.cpp src/app/BondAngleSuggester.h src/app/BondAngleSuggester.cpp)
target_include_directories(bond_angle_suggester_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
)
target_link_libraries(bond_angle_suggester_test PRIVATE Qt6::Core)
add_test(NAME bond_angle_suggester_test COMMAND bond_angle_suggester_test)

qt_add_executable(document_state_test tests/document_state_test.cpp)
target_link_libraries(document_state_test PRIVATE molecule_model)
target_compile_definitions(document_state_test PRIVATE SKETCH_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME document_state_test COMMAND document_state_test)

qt_add_executable(render_primitives_test tests/render_primitives_test.cpp)
target_link_libraries(render_primitives_test PRIVATE molecule_model)
target_compile_definitions(render_primitives_test PRIVATE SKETCH_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME render_primitives_test COMMAND render_primitives_test)

qt_add_executable(sdf_batch_test tests/sdf_batch_test.cpp)
target_link_libraries(sdf_batch_test PRIVATE molecule_model)
target_compile_definitions(sdf_batch_test PRIVATE SKETCH_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME sdf_batch_test COMMAND sdf_batch_test)

qt_add_executable(clipboard_preview_test tests/clipboard_preview_test.cpp)
target_link_libraries(clipboard_preview_test PRIVATE molecule_model)
target_compile_definitions(clipboard_preview_test PRIVATE SKETCH_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME clipboard_preview_test COMMAND clipboard_preview_test)

qt_add_executable(biopolymer_sequence_view_test tests/biopolymer_sequence_view_test.cpp)
target_link_libraries(biopolymer_sequence_view_test PRIVATE molecule_model)
add_test(NAME biopolymer_sequence_view_test COMMAND biopolymer_sequence_view_test)

qt_add_executable(render_primitives_to_variant_test tests/render_primitives_to_variant_test.cpp)
target_link_libraries(render_primitives_to_variant_test PRIVATE molecule_model)
add_test(NAME render_primitives_to_variant_test COMMAND render_primitives_to_variant_test)

