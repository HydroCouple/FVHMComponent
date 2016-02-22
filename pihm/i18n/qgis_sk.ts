<!DOCTYPE TS><TS>
<context>
    <name>CoordinateCapture</name>
    <message>
        <source>Coordinate Capture</source>
        <translation>Odchyt súradníc</translation>
    </message>
    <message>
        <source>Click on the map to view coordinates and capture to clipboard.</source>
        <translation>Kliknutím na mapu sa zobrazia súradnice bodu kliknutia a a tieto zároveň budú skopírované do schránky.</translation>
    </message>
    <message>
        <source>&amp;Coordinate Capture</source>
        <translation>&amp;Odchyt súradníc</translation>
    </message>
    <message>
        <source>Copy to clipboard</source>
        <translation>Kopírovať do schránky</translation>
    </message>
    <message>
        <source>Click to select the CRS to use for coordinate display</source>
        <translation>Kliknutím možno vybrať súradnicový systém, v ktorom budú súradnice zobrazené</translation>
    </message>
    <message>
        <source>Coordinate in your selected CRS</source>
        <translation>Súradnice vo vybranom CRS</translation>
    </message>
    <message>
        <source>Coordinate in map canvas coordinate reference system</source>
        <translation>Súradnice v súradnicovom systéme mapového plátna</translation>
    </message>
    <message>
        <source>Click to enable mouse tracking. Click the canvas to stop</source>
        <translation>Kliknutím sa zapne sledovanie myši. Kliknutím na plátno bude sledovanie vypnuté</translation>
    </message>
</context>
<context>
    <name>CoordinateCaptureGui</name>
    <message>
        <source>Welcome to your automatically generated plugin!</source>
        <translation>Vitajte vo vašom automaticky vygenerovanom zásuvnom module!</translation>
    </message>
    <message>
        <source>This is just a starting point. You now need to modify the code to make it do something useful....read on for a more information to get yourself started.</source>
        <translation>Toto je len začiatok. Teraz je potrebné upraviť kód, tak aby robil niečo užitočné...čítajte ďalej a dozviete sa viac informácií, ktoré Vám pomôžu rozbehnúť sa.</translation>
    </message>
    <message>
        <source>Documentation:</source>
        <translation>Dokumentácia:</translation>
    </message>
    <message>
        <source>You really need to read the QGIS API Documentation now at:</source>
        <translation type="unfinished">Teraz je potrebné si prečítať dokumentáciu k API rozhraniu aplikácie QGIS na:</translation>
    </message>
    <message>
        <source>In particular look at the following classes:</source>
        <translation type="unfinished">Zvlášť sa pozrite na nasledujúce triedy:</translation>
    </message>
    <message>
        <source>QgsPlugin is an ABC that defines required behaviour your plugin must provide. See below for more details.</source>
        <translation type="unfinished">QgsPlugin je ABC, ktorá definuje požadované správanie sa zásuvného modulu, ktoré musí zásvuný modul zabezpečovať. Viac informácií sa dočítate ďalej.</translation>
    </message>
    <message>
        <source>What are all the files in my generated plugin directory for?</source>
        <translation type="unfinished">Na čo slúžia všetky tie súbory v adresári vygenerovanom pre môj zásuvný modul?</translation>
    </message>
    <message>
        <source>This is the generated CMake file that builds the plugin. You should add you application specific dependencies and source files to this file.</source>
        <translation type="unfinished">To je vygenerovný súbor CMake ktorý slúži na kompiláciu zásuvného modulu. Do neho sa pridávajú závislosti pre vašu aplikáciu a zdrojové súbory.</translation>
    </message>
    <message>
        <source>This is the class that provides the &apos;glue&apos; between your custom application logic and the QGIS application. You will see that a number of methods are already implemented for you - including some examples of how to add a raster or vector layer to the main application map canvas. This class is a concrete instance of the QgisPlugin interface which defines required behaviour for a plugin. In particular, a plugin has a number of static methods and members so that the QgsPluginManager and plugin loader logic can identify each plugin, create an appropriate menu entry for it etc. Note there is nothing stopping you creating multiple toolbar icons and menu entries for a single plugin. By default though a single menu entry and toolbar button is created and its pre-configured to call the run() method in this class when selected. This default implementation provided for you by the plugin builder is well documented, so please refer to the code for further advice.</source>
        <translation type="unfinished">Toto je trieda ktorá zabezpečuje rozhranie medzi vlastnou logikou aplikácie a aplikáciou QGIS. Množstvo metód je už implementovaných pre vás vrátane niektorých príkladov ako pridať rastrovú alebo vektorovú vrstvu na mapové plátno hlavnej aplikácie. Táto trieda je konkrétnou inštanciou rozhrania QgisPlugin ktoré definuje správanie zásuvného modulu. Napríklad zásuvný modul má množstvo statických metód a členov takže QgsPluginManager a logika nahrávania zásuvných modulov môže identifikovať každý zásuvný modul, vytvoriť preň príslušné položky menu atď. Nie je nič, čo by vám bránilo vo vytváraní viacerých ikôn na panely a položiek  v menu pre jeden zásuvný modul. Štandardne je vytvorená jedna položka menu a tlačidlo na panely a je prednastavená na volanie metódy run() v tejto triede keď je vybratá. Táto predvolená implementácia pripravená pre vás builderom zásuvných modulov je dobre dokumentovaná takže ďalšie rady získate priamo v kóde.</translation>
    </message>
    <message>
        <source>This is a Qt designer &apos;ui&apos; file. It defines the look of the default plugin dialog without implementing any application logic. You can modify this form to suite your needs or completely remove it if your plugin does not need to display a user form (e.g. for custom MapTools).</source>
        <translation type="unfinished">Toto je &apos;ui&apos; súbor Qt dizajnéra. Definuje vzhľad štandardného dialógového okna zásuvného modulu bez implementácie akejkoľvek aplikačnej logiky. Tento formulár môžete upraviť tak, aby vyhovoval vašim potrebám, alebo ho úplne odobrať ak váš zásuvný modul nepotrebuje zobraziť formulár pre používateľa (napr. pre vlastné MapTools).</translation>
    </message>
    <message>
        <source>This is the concrete class where application logic for the above mentioned dialog should go. The world is your oyster here really....</source>
        <translation type="unfinished">Toto je konkrétna trieda, kde by mala byť aplikačná logika pre vyššie spomenuté dialógové okno. Svet je váš ?? tu skutočne...</translation>
    </message>
    <message>
        <source>This is the Qt4 resources file for your plugin. The Makefile generated for your plugin is all set up to compile the resource file so all you need to do is add your additional icons etc using the simple xml file format. Note the namespace used for all your resources e.g. (&apos;:/Homann/&apos;). It is important to use this prefix for all your resources. We suggest you include any other images and run time data in this resurce file too.</source>
        <translation type="unfinished">Toto je súbor zdrojov Qt4 pre váš zásuvný modul. Makefile generovaný pre váš zásuvný modul je nastavený aby skompiloval súbor zdrojov, takže všetko čo potrebujete spraviť je pridať vaše ďalšie ikony atď. použitím jednoduchého súboru vo formáte xml. Nezabudnite na menný priestor použitý pre všetky vaše zdroje napr (&apos;:/Homann/&apos;). Je dôležité použiť predponu pre všetky vaše zdroje. Odporúča sa zahrnúť všetky ďalšie obrázky a run time údaje do tohoto súboru zdrojov tiež.</translation>
    </message>
    <message>
        <source>This is the icon that will be used for your plugin menu entry and toolbar icon. Simply replace this icon with your own icon to make your plugin disctinctive from the rest.</source>
        <translation type="unfinished">Toto je ikona pre váš zásuvný modul, ktorá bude použitá pre menu a ikonu na paneli. Jednoducho zameňte túto ikonu s vašou vlastnou ikonou aby bol váš zásuvný modul odlíšiteľný od ostatných.</translation>
    </message>
    <message>
        <source>This file contains the documentation you are reading now!</source>
        <translation type="unfinished">Tento súbor obsahuje dokumentáciu, ktorú práve čítate!</translation>
    </message>
    <message>
        <source>Getting developer help:</source>
        <translation>Ako získať pomoc od vývojárov:</translation>
    </message>
    <message>
        <source>For Questions and Comments regarding the plugin builder template and creating your features in QGIS using the plugin interface please contact us via:</source>
        <translation type="unfinished">Pre otázky a komentáre súvisiace so šablónou zásuvného modulu a vytváraním vlastnej funkcionality v QGISe s použitím rozhrania zásuvného modulu nás kontaktujte cez:</translation>
    </message>
    <message>
        <source>&lt;li&gt; the QGIS developers mailing list, or &lt;/li&gt;&lt;li&gt; IRC (#qgis on freenode.net)&lt;/li&gt;</source>
        <translation type="unfinished">&lt;li&gt; e-mailová konferencia vývojárov QGISu, alebo &lt;/li&gt;&lt;li&gt; IRC (kanál #qgis na freenode.net)&lt;/li&gt;</translation>
    </message>
    <message>
        <source>QGIS is distributed under the Gnu Public License. If you create a useful plugin please consider contributing it back to the community.</source>
        <translation type="unfinished">QGIS je šírený pod licenciou Gnu Public License. Ak vytvoríte užitočný zásuvný modul, prosím prispejte s ním do komunity.</translation>
    </message>
    <message>
        <source>Have fun and thank you for choosing QGIS.</source>
        <translation type="unfinished">Prajeme Vám veľa zábavy s aplikáciou QGIS.</translation>
    </message>
</context>
<context>
    <name>CoordinateCaptureGuiBase</name>
    <message>
        <source>QGIS Plugin Template</source>
        <translation>QGIS šablóna zásuvného modulu</translation>
    </message>
    <message>
        <source>Plugin Template</source>
        <translation>Šablóna zásuvného modulu</translation>
    </message>
</context>
<context>
    <name>Dialog</name>
    <message>
        <source>QGIS Plugin Installer</source>
        <translation type="obsolete">QGIS Inštalátor zásuvných modulov</translation>
    </message>
    <message>
        <source>Name of plugin to install</source>
        <translation type="obsolete">Bude nainštalovaný zásuvný modul:</translation>
    </message>
    <message>
        <source>Get List</source>
        <translation type="obsolete">Získať zoznam</translation>
    </message>
    <message>
        <source>Done</source>
        <translation type="obsolete">Dokončiť</translation>
    </message>
    <message>
        <source>Install Plugin</source>
        <translation type="obsolete">Inštalovať zásuvný modul</translation>
    </message>
    <message>
        <source>The plugin will be installed to ~/.qgis/python/plugins</source>
        <translation type="obsolete">Zásuvný modul bude nainštalovaný do ~/.qgis/python/plugins</translation>
    </message>
    <message>
        <source>Name</source>
        <translation type="obsolete">Meno</translation>
    </message>
    <message>
        <source>Version</source>
        <translation type="obsolete">Verzia</translation>
    </message>
    <message>
        <source>Description</source>
        <translation type="obsolete">Popis</translation>
    </message>
    <message>
        <source>Author</source>
        <translation type="obsolete">Autor</translation>
    </message>
    <message>
        <source>Select repository, retrieve the list of available plugins, select one and install it</source>
        <translation type="obsolete">Vyberte repozitár, získajte zoznam dostupných zásuvných modulov, vyberte jeden a nainštalujte ho</translation>
    </message>
    <message>
        <source>Repository</source>
        <translation type="obsolete">Repozitár</translation>
    </message>
    <message>
        <source>Active repository:</source>
        <translation type="obsolete">Aktívny repozitár:</translation>
    </message>
    <message>
        <source>Add</source>
        <translation type="obsolete">Pridať</translation>
    </message>
    <message>
        <source>Edit</source>
        <translation type="obsolete">Upraviť</translation>
    </message>
    <message>
        <source>Delete</source>
        <translation type="obsolete">Vymazať</translation>
    </message>
    <message>
        <source>Connect</source>
        <translation>Spojiť</translation>
    </message>
    <message>
        <source>Browse</source>
        <translation>Prechádzať</translation>
    </message>
    <message>
        <source>OGR Converter</source>
        <translation>OGR Prevodník</translation>
    </message>
    <message>
        <source>Could not establish connection to: &apos;</source>
        <translation>Nemožno nadviazať spojenie so serverom: &apos;</translation>
    </message>
    <message>
        <source>Open OGR file</source>
        <translation>Otvoriť súbor podporovaný knižnicou OGR</translation>
    </message>
    <message>
        <source>OGR File Data Source (*.*)</source>
        <translation type="unfinished">Súbory OGR zdroje údajov (*.*)</translation>
    </message>
    <message>
        <source>Open Directory</source>
        <translation>Otvoriť adresár</translation>
    </message>
    <message>
        <source>Input OGR dataset is missing!</source>
        <translation type="unfinished">Chýba vstupná súprava údajov OGR!</translation>
    </message>
    <message>
        <source>Input OGR layer name is missing!</source>
        <translation>Chýba meno vstupnej OGR vrstvy!</translation>
    </message>
    <message>
        <source>Target OGR format not selected!</source>
        <translation>Nie je vybraný cieľový formát OGR !</translation>
    </message>
    <message>
        <source>Output OGR dataset is missing!</source>
        <translation type="unfinished">Chýba výstupná súprava údajov OGR!</translation>
    </message>
    <message>
        <source>Output OGR layer name is missing!</source>
        <translation>Chýba meno výstupnej OGR vrstvy!</translation>
    </message>
    <message>
        <source>Successfully translated layer &apos;</source>
        <translation>Úspešne preložená vrstva &apos;</translation>
    </message>
    <message>
        <source>Failed to translate layer &apos;</source>
        <translation>Zlyhalo preloženie vrstvy &apos;</translation>
    </message>
    <message>
        <source>Successfully connected to: &apos;</source>
        <translation>Úspešne pripojené k: &apos;</translation>
    </message>
    <message>
        <source>Choose a file name to save to</source>
        <translation type="unfinished">Vyberte meno súboru do ktorého sa bude ukladať</translation>
    </message>
</context>
<context>
    <name>Gui</name>
    <message>
        <source>Welcome to your automatically generated plugin!</source>
        <translation>Vitajte v automaticky vygenerovanom zásuvnom module!</translation>
    </message>
    <message>
        <source>This is just a starting point. You now need to modify the code to make it do something useful....read on for a more information to get yourself started.</source>
        <translation>Toto je len začiatok. Teraz je potrebné upraviť kód, tak aby robil niečo užitočné...čítajte ďalej a dozviete sa viac informácií, ktoré Vám pomôžu rozbehnúť sa.</translation>
    </message>
    <message>
        <source>Documentation:</source>
        <translation>Dokumentácia:</translation>
    </message>
    <message>
        <source>You really need to read the QGIS API Documentation now at:</source>
        <translation>Teraz je potrebné si prečítať dokumentáciu k API rozhraniu aplikácie QGIS na:</translation>
    </message>
    <message>
        <source>In particular look at the following classes:</source>
        <translation type="unfinished">Zvlášť sa pozrite na nasledujúce triedy:</translation>
    </message>
    <message>
        <source>QgsPlugin is an ABC that defines required behaviour your plugin must provide. See below for more details.</source>
        <translation type="unfinished">QgsPlugin je ABC, ktorá definuje požadované správanie sa zásuvného modulu, ktoré musí zabezpečiť. Viac informácií sa dočítate ďalej.</translation>
    </message>
    <message>
        <source>What are all the files in my generated plugin directory for?</source>
        <translation type="unfinished">Na čo slúžia všetky tie súbory v adresári vygenerovanom pre môj zásuvný modul?</translation>
    </message>
    <message>
        <source>This is the generated CMake file that builds the plugin. You should add you application specific dependencies and source files to this file.</source>
        <translation type="unfinished">To je vygenerovný súbor CMake ktorý slúži na kompiláciu zásuvného modulu. Do neho sa pridávajú závislosti pre vašu aplikáciu a zdrojové súbory.</translation>
    </message>
    <message>
        <source>This is the class that provides the &apos;glue&apos; between your custom application logic and the QGIS application. You will see that a number of methods are already implemented for you - including some examples of how to add a raster or vector layer to the main application map canvas. This class is a concrete instance of the QgisPlugin interface which defines required behaviour for a plugin. In particular, a plugin has a number of static methods and members so that the QgsPluginManager and plugin loader logic can identify each plugin, create an appropriate menu entry for it etc. Note there is nothing stopping you creating multiple toolbar icons and menu entries for a single plugin. By default though a single menu entry and toolbar button is created and its pre-configured to call the run() method in this class when selected. This default implementation provided for you by the plugin builder is well documented, so please refer to the code for further advice.</source>
        <translation type="unfinished">Toto je trieda ktorá zabezpečuje rozhranie medzi vlastnou logikou aplikácie a aplikáciou QGIS. Množstvo metód je už implementovaných pre vás vrátane niektorých príkladov ako pridať rastrovú alebo vektorovú vrstvu na mapové plátno hlavnej aplikácie. Táto trieda je konkrétnou inštanciou rozhrania QgisPlugin ktoré definuje správanie zásuvného modulu. Napríklad zásuvný modul má množstvo statických metód a členov takže QgsPluginManager a logika nahrávania zásuvných modulov môže identifikovať každý zásuvný modul, vytvoriť preň príslušné položky menu atď. Nie je nič, čo by vám bránilo vo vytváraní viacerých ikôn na panely a položiek  v menu pre jeden zásuvný modul. Štandardne je vytvorená jedna položka menu a tlačidlo na panely a je prednastavená na volanie metódy run() v tejto triede keď je vybratá. Táto predvolená implementácia pripravená pre vás builderom zásuvných modulov je dobre dokumentovaná takže ďalšie rady získate priamo v kóde.</translation>
    </message>
    <message>
        <source>This is a Qt designer &apos;ui&apos; file. It defines the look of the default plugin dialog without implementing any application logic. You can modify this form to suite your needs or completely remove it if your plugin does not need to display a user form (e.g. for custom MapTools).</source>
        <translation type="unfinished">Toto je &apos;ui&apos; súbor Qt dizajnéra. Definuje vzhľad štandardného dialógového okna zásuvného modulu bez implementácie akejkoľvek aplikačnej logiky. Tento formulár môžete tak aby vyhovoval vašim potrebám, alebo ho úplne odobrať ak váš zásuvný modul nepotrebuje zobraziť formulár pre používateľa (napr. pre vlastné MapTools).</translation>
    </message>
    <message>
        <source>This is the concrete class where application logic for the above mentioned dialog should go. The world is your oyster here really....</source>
        <translation type="unfinished">Toto je konkrétna trieda, kde by mala byť aplikačná logika pre vyššie spomenuté dialógové okno. Svet je váš ?? tu skutočne...</translation>
    </message>
    <message>
        <source>This is the Qt4 resources file for your plugin. The Makefile generated for your plugin is all set up to compile the resource file so all you need to do is add your additional icons etc using the simple xml file format. Note the namespace used for all your resources e.g. (&apos;:/Homann/&apos;). It is important to use this prefix for all your resources. We suggest you include any other images and run time data in this resurce file too.</source>
        <translation type="unfinished">Toto je súbor zdrojov Qt4 pre váš zásuvný modul. Makefile generovaný pre váš zásuvný modul je nastavený aby skompiloval súbor zdrojov, takže všetko čo potrebujete spraviť je pridať vaše ďalšie ikony atď. použitím jednoduchého súboru vo formáte xml. Nezabudnite na menný priestor použitý pre všetky vaše zdroje napr (&apos;:/Homann/&apos;). Je dôležité použiť predponu pre všetky vaše zdroje. Odporúča sa zahrnúť všetky ďalšie obrázky a run time údaje do tohoto súboru zdrojov tiež.</translation>
    </message>
    <message>
        <source>This is the icon that will be used for your plugin menu entry and toolbar icon. Simply replace this icon with your own icon to make your plugin disctinctive from the rest.</source>
        <translation type="unfinished">Toto je ikona pre váš zásuvný modul, ktorá bude použitá pre menu a ikonu na paneli. Jednoducho zameňte túto ikonu s vašou vlastnou ikonou aby bol váš zásuvný modul odlíšiteľný od ostatných.</translation>
    </message>
    <message>
        <source>This file contains the documentation you are reading now!</source>
        <translation>Tento súbor obsahuje dokumentáciu, ktorú práve čítate!</translation>
    </message>
    <message>
        <source>Getting developer help:</source>
        <translation>Získanie pomoci od vývojárov:</translation>
    </message>
    <message>
        <source>For Questions and Comments regarding the plugin builder template and creating your features in QGIS using the plugin interface please contact us via:</source>
        <translation type="unfinished">Pre otázky a komentáre súvisiace so šablónou záusvného modulu a vytváraním vlastnej funkcionality v QGISe s použitím rozhrania zásuvného modulu nás kontaktujte cez:</translation>
    </message>
    <message>
        <source>&lt;li&gt; the QGIS developers mailing list, or &lt;/li&gt;&lt;li&gt; IRC (#qgis on freenode.net)&lt;/li&gt;</source>
        <translation>&lt;li&gt; e-mailová konferencia vývojárov QGISu, alebo &lt;/li&gt;&lt;li&gt; IRC (kanál #qgis na freenode.net)&lt;/li&gt;</translation>
    </message>
    <message>
        <source>QGIS is distributed under the Gnu Public License. If you create a useful plugin please consider contributing it back to the community.</source>
        <translation>QGIS je šírený pod licenciou Gnu Public License. Ak vytvoríte užitočný zásuvný modul, prosím prispejte ním do komunity.</translation>
    </message>
    <message>
        <source>Have fun and thank you for choosing QGIS.</source>
        <translation>Ďakujeme, že ste si vybrali QGIS a prajeme veľa zábavy.</translation>
    </message>
</context>
<context>
    <name>MapCoordsDialogBase</name>
    <message>
        <source>Enter map coordinates</source>
        <translation>Zadanie mapových súradníc</translation>
    </message>
    <message>
        <source>X:</source>
        <translation>X:</translation>
    </message>
    <message>
        <source>Y:</source>
        <translation>Y:</translation>
    </message>
    <message>
        <source>&amp;OK</source>
        <translation>&amp;OK</translation>
    </message>
    <message>
        <source>&amp;Cancel</source>
        <translation>&amp;Zrušiť</translation>
    </message>
    <message>
        <source>Enter X and Y coordinates which correspond with the selected point on the image. Alternatively, click the button with icon of a pencil and then click a corresponding point on map canvas of QGIS to fill in coordinates of that point.</source>
        <translation>Zadajte X-ovú a Y-ovú súradnicu patriacu vybranému bodu na snímke. Druhou možnosťou ako ich zadať je kliknúť na ikonu s ceruzkou a následným kliknutím na zodpovedajúci bod na mapovom plátne QGISu sa vyplnia políčka súradnice daného bodu.</translation>
    </message>
    <message>
        <source> from map canvas</source>
        <translation> Z mapového plátna</translation>
    </message>
</context>
<context>
    <name>OgrConverterGuiBase</name>
    <message>
        <source>OGR Layer Converter</source>
        <translation>Prevodník vrstiev OGR</translation>
    </message>
    <message>
        <source>Source</source>
        <translation>Zdroj</translation>
    </message>
    <message>
        <source>Format</source>
        <translation>Formát</translation>
    </message>
    <message>
        <source>File</source>
        <translation>Súbor</translation>
    </message>
    <message>
        <source>Directory</source>
        <translation>Adresár</translation>
    </message>
    <message>
        <source>Remote source</source>
        <translation>Vzdialený zdroj</translation>
    </message>
    <message>
        <source>Dataset</source>
        <translation type="unfinished">Súprava údajov</translation>
    </message>
    <message>
        <source>Browse</source>
        <translation>Prechádzať</translation>
    </message>
    <message>
        <source>Layer</source>
        <translation>Vrstva</translation>
    </message>
    <message>
        <source>Target</source>
        <translation>Cieľ</translation>
    </message>
</context>
<context>
    <name>OgrPlugin</name>
    <message>
        <source>Run OGR Layer Converter</source>
        <translation>Spustiť prevodník OGR vrstiev</translation>
    </message>
    <message>
        <source>Translates vector layers between formats supported by OGR library</source>
        <translation>Prevádza vektorové vrstvy medzi rôznymi formátmi podporovanými knižnicou OGR</translation>
    </message>
    <message>
        <source>OG&amp;R Converter</source>
        <translation>Prevodník OG&amp;R </translation>
    </message>
</context>
<context>
    <name>QFileDialog</name>
    <message>
        <source>Save experiment report to portable document format (.pdf)</source>
        <translation type="unfinished">Uložiť správu z experimentu do formátu PDF (.pdf)</translation>
    </message>
    <message>
        <source>Load layer properties from style file (.qml)</source>
        <translation>Nahrať vlastnosti vrstvy zo súboru štýlu (.qml)</translation>
    </message>
    <message>
        <source>Save layer properties as style file (.qml)</source>
        <translation>Uložiť vlastnosti vrstvy do súboru štýlu (.qml)</translation>
    </message>
</context>
<context>
    <name>QObject</name>
    <message>
        <source>No Data Providers</source>
        <translation>Žiadne nástroje na prístup k údajom</translation>
    </message>
    <message>
        <source>No Data Provider Plugins</source>
        <comment>No QGIS data provider plugins found in:</comment>
        <translation>Žiadne zásuvné moduly na prístup k údajom</translation>
    </message>
    <message>
        <source>No vector layers can be loaded. Check your QGIS installation</source>
        <translation>Nie je možné nahrať žiadne vektorové vrstvy. Skontrolujte vašu inštaláciu QGIS</translation>
    </message>
    <message>
        <source>No data provider plugins are available. No vector layers can be loaded</source>
        <translation>Žiadne zásuvné moduly na prístup k údajom nie sú k dispozícii. Žiadne vektorové vrstvy nemôžu byť nahraté</translation>
    </message>
    <message>
        <source>QGis files (*.qgs)</source>
        <translation>QGIS súbory (*.qgs)</translation>
    </message>
    <message>
        <source> at line </source>
        <translation> v riadku </translation>
    </message>
    <message>
        <source> column </source>
        <translation> stĺpec </translation>
    </message>
    <message>
        <source> for file </source>
        <translation> pre súbor </translation>
    </message>
    <message>
        <source>Unable to save to file </source>
        <translation>Nie je možné uložiť súbor</translation>
    </message>
    <message>
        <source>Referenced column wasn&apos;t found: </source>
        <translation> Stĺpec na ktorý smeroval odkaz sa nenašiel: </translation>
    </message>
    <message>
        <source>Division by zero.</source>
        <translation>Delenie nulou.</translation>
    </message>
    <message>
        <source>No active layer</source>
        <translation>Žiadna vrstva nie je aktívna</translation>
    </message>
    <message>
        <source>Band</source>
        <translation>Kanál</translation>
    </message>
    <message>
        <source>action</source>
        <translation>akcia</translation>
    </message>
    <message>
        <source> features found</source>
        <translation> objektov nájdených</translation>
    </message>
    <message>
        <source> 1 feature found</source>
        <translation> 1 objekt nájdený</translation>
    </message>
    <message>
        <source>No features found</source>
        <translation>Nenašli sa žiadne objekty</translation>
    </message>
    <message>
        <source>No features were found in the active layer at the point you clicked</source>
        <translation>V bode, na ktorý ste klikli sa v aktívnej vrstve nenašli žiadne objekty </translation>
    </message>
    <message>
        <source>Could not identify objects on</source>
        <translation>Nemožno identifikovať objekty na</translation>
    </message>
    <message>
        <source>because</source>
        <translation>z nasledujúceho dôvodu</translation>
    </message>
    <message>
        <source>New centroid</source>
        <translation>Nový centroid</translation>
    </message>
    <message>
        <source>New point</source>
        <translation>Nový bod</translation>
    </message>
    <message>
        <source>New vertex</source>
        <translation>Nový uzol</translation>
    </message>
    <message>
        <source>Undo last point</source>
        <translation>Posledný bod späť</translation>
    </message>
    <message>
        <source>Close line</source>
        <translation>Uzavrieť líniu</translation>
    </message>
    <message>
        <source>Select vertex</source>
        <translation>Vybrať uzol</translation>
    </message>
    <message>
        <source>Select new position</source>
        <translation>Vybrať novú polohu</translation>
    </message>
    <message>
        <source>Select line segment</source>
        <translation>Vybrať líniový úsek</translation>
    </message>
    <message>
        <source>New vertex position</source>
        <translation>Nová poloha uzla</translation>
    </message>
    <message>
        <source>Release</source>
        <translation>Uvoľniť</translation>
    </message>
    <message>
        <source>Delete vertex</source>
        <translation>Vymazať uzol</translation>
    </message>
    <message>
        <source>Release vertex</source>
        <translation>Uvoľniť uzol</translation>
    </message>
    <message>
        <source>Select element</source>
        <translation>Vybrať element</translation>
    </message>
    <message>
        <source>New location</source>
        <translation>Nová lokalita (location)</translation>
    </message>
    <message>
        <source>Release selected</source>
        <translation>Uvoľniť vybrané</translation>
    </message>
    <message>
        <source>Delete selected / select next</source>
        <translation>Vymazať vybrané / vybrať ďalšie</translation>
    </message>
    <message>
        <source>Select position on line</source>
        <translation>Vybrať polohu na línii</translation>
    </message>
    <message>
        <source>Split the line</source>
        <translation>Rozdeliť líniu</translation>
    </message>
    <message>
        <source>Release the line</source>
        <translation>Uvoľniť líniu</translation>
    </message>
    <message>
        <source>Select point on line</source>
        <translation>Vybrať bod na línii</translation>
    </message>
    <message>
        <source>Length</source>
        <translation>Dĺžka</translation>
    </message>
    <message>
        <source>Area</source>
        <translation>Rozloha</translation>
    </message>
    <message>
        <source>Label</source>
        <translation>Popis</translation>
    </message>
    <message>
        <source>Project file read error: </source>
        <translation> Chyba pri čítaní súboru projketu:</translation>
    </message>
    <message>
        <source>Fit to a linear transform requires at least 2 points.</source>
        <translation> Na lineárnu transformáciu sú potrebné najmenej dva body.</translation>
    </message>
    <message>
        <source>Fit to a Helmert transform requires at least 2 points.</source>
        <translation>Na Helmertovu transformáciu sú potrebné najmenej dva body.</translation>
    </message>
    <message>
        <source>Fit to an affine transform requires at least 4 points.</source>
        <translation>Na afinnú transformáciu sú potrebné najmenej štyri body.</translation>
    </message>
    <message>
        <source>Couldn&apos;t open the data source: </source>
        <translation> Nemožno otvoriť zdroj údajov:</translation>
    </message>
    <message>
        <source>Parse error at line </source>
        <translation> Chyba pri parse v riadku</translation>
    </message>
    <message>
        <source>GPS eXchange format provider</source>
        <translation>Správca údajov vo formáte GPS eXchange</translation>
    </message>
    <message>
        <source>Caught a coordinate system exception while trying to transform a point. Unable to calculate line length.</source>
        <translation>Pri pokuse o transformáciu bodu bola zachytená výnimka súradnicového systému. Dĺžku línie nemožno vypočítať.</translation>
    </message>
    <message>
        <source>Caught a coordinate system exception while trying to transform a point. Unable to calculate polygon area.</source>
        <translation>Pri pokuse o transformáciu bodu bola zachytená výnimka súradnicového systému. Plochu polygónu nemožno vypočítať.</translation>
    </message>
    <message>
        <source>CopyrightLabel</source>
        <translation>Označenie autorských práv</translation>
    </message>
    <message>
        <source>Draws copyright information</source>
        <translation>Vypíše na plátno informáciu o autorských právach</translation>
    </message>
    <message>
        <source>Version 0.1</source>
        <translation>Verzia 0.1</translation>
    </message>
    <message>
        <source>Version 0.2</source>
        <translation>Verzia 0.2</translation>
    </message>
    <message>
        <source>Loads and displays delimited text files containing x,y coordinates</source>
        <translation>Nahrá a zobrazí súbory s oddeleným textom obsahujúce súradnice x a y</translation>
    </message>
    <message>
        <source>Add Delimited Text Layer</source>
        <translation>Pridanie vrstvy z oddeleného textu</translation>
    </message>
    <message>
        <source>Georeferencer</source>
        <translation>Georeferencer</translation>
    </message>
    <message>
        <source>Adding projection info to rasters</source>
        <translation>Pridáva rastrom informácie o mapovom zobrazení</translation>
    </message>
    <message>
        <source>GPS Tools</source>
        <translation>Nástroje na prácu s GPS</translation>
    </message>
    <message>
        <source>Tools for loading and importing GPS data</source>
        <translation>Nástroje na nahratie a import údajov GPS</translation>
    </message>
    <message>
        <source>GRASS</source>
        <translation>GRASS</translation>
    </message>
    <message>
        <source>GRASS layer</source>
        <translation>vrstva GRASSu</translation>
    </message>
    <message>
        <source>Graticule Creator</source>
        <translation>Tvorba súradnicovej siete</translation>
    </message>
    <message>
        <source>Builds a graticule</source>
        <translation>Vytvorí súradnicovú sieť</translation>
    </message>
    <message>
        <source>NorthArrow</source>
        <translation>Smerová ružica</translation>
    </message>
    <message>
        <source>Displays a north arrow overlayed onto the map</source>
        <translation>Zobrazí na mape smerovú ružicu</translation>
    </message>
    <message>
        <source>[menuitemname]</source>
        <translation>[názovpoložkymenu]</translation>
    </message>
    <message>
        <source>[plugindescription]</source>
        <translation>[popiszásuvnéhomodulu]</translation>
    </message>
    <message>
        <source>ScaleBar</source>
        <translation>Grafická mierka</translation>
    </message>
    <message>
        <source>Draws a scale bar</source>
        <translation>Vykreslí grafickú mierku</translation>
    </message>
    <message>
        <source>SPIT</source>
        <translation>SPIT</translation>
    </message>
    <message>
        <source>Shapefile to PostgreSQL/PostGIS Import Tool</source>
        <translation>Nástroj na import súborov Shape do PostgreSQL/PostGIS</translation>
    </message>
    <message>
        <source>WFS plugin</source>
        <translation>Zásuvný modul WFS</translation>
    </message>
    <message>
        <source>Adds WFS layers to the QGIS canvas</source>
        <translation>Do mapy pridá vrstvy WFS</translation>
    </message>
    <message>
        <source>GRASS plugin</source>
        <translation>Zásuvný modul GRASSu</translation>
    </message>
    <message>
        <source>QGIS couldn&apos;t find your GRASS installation.
Would you like to specify path (GISBASE) to your GRASS installation?</source>
        <translation>QGIS nemôže nájsť vašu inštaláciu GRASSu.
Chcete zadať cestu (GISBASE) k vašej inštalácii GRASSu?</translation>
    </message>
    <message>
        <source>Choose GRASS installation path (GISBASE)</source>
        <translation type="unfinished">Vyberte cestu k inštalácii GRASSu (GISBASE)</translation>
    </message>
    <message>
        <source>GRASS data won&apos;t be available if GISBASE is not specified.</source>
        <translation>Údaje GRASSu nebudú dostupné pokiaľ nebude určená GISBASE.</translation>
    </message>
    <message>
        <source>Not a vector layer</source>
        <translation>Nie je vektorová vrstva</translation>
    </message>
    <message>
        <source>The current layer is not a vector layer</source>
        <translation>Táto vrstva nie je vektorovou vrstvou</translation>
    </message>
    <message>
        <source>Layer not editable</source>
        <translation>Vrstva nie je upravovateľná</translation>
    </message>
    <message>
        <source>Cannot edit the vector layer. To make it editable, go to the file item of the layer, right click and check &apos;Allow Editing&apos;.</source>
        <translation>Túto vektorovú vrstvu nemožno upravovať. Aby ju bolo možné upravovať je potrebné kliknúť pravý tlačidlom myši na meno vrstvy a zaškrtnúť &apos;Povoliť úpravy&apos;.</translation>
    </message>
    <message>
        <source>Wrong editing tool</source>
        <translation>Nesprávny naástroj na úpravy</translation>
    </message>
    <message>
        <source>Cannot apply the &apos;capture point&apos; tool on this vector layer</source>
        <translation>Na túto vektorovú vrstvu nemožno používať nástroj &apos;Ziskať bod&apos;</translation>
    </message>
    <message>
        <source>Cannot apply the &apos;capture line&apos; tool on this vector layer</source>
        <translation>Na túto vektorovú vrstvu nemožno používať nástroj &apos;Ziskať líniu&apos;</translation>
    </message>
    <message>
        <source>Cannot apply the &apos;capture polygon&apos; tool on this vector layer</source>
        <translation>Na túto vektorovú vrstvu nemožno používať nástroj &apos;Ziskať polygón&apos;</translation>
    </message>
    <message>
        <source>Error</source>
        <translation>Chyba</translation>
    </message>
    <message>
        <source> km2</source>
        <translation type="unfinished"> km2</translation>
    </message>
    <message>
        <source> ha</source>
        <translation> ha</translation>
    </message>
    <message>
        <source> m2</source>
        <translation type="unfinished"> m2</translation>
    </message>
    <message>
        <source> m</source>
        <translation> m</translation>
    </message>
    <message>
        <source> km</source>
        <translation> km</translation>
    </message>
    <message>
        <source> mm</source>
        <translation> mm</translation>
    </message>
    <message>
        <source> cm</source>
        <translation> cm</translation>
    </message>
    <message>
        <source> sq mile</source>
        <translation> štv. míľ</translation>
    </message>
    <message>
        <source> sq ft</source>
        <translation> štv. stôp</translation>
    </message>
    <message>
        <source> mile</source>
        <translation> míľ</translation>
    </message>
    <message>
        <source> foot</source>
        <translation> stopa</translation>
    </message>
    <message>
        <source> feet</source>
        <translation> stopy</translation>
    </message>
    <message>
        <source> sq.deg.</source>
        <translation> štv. stupň.</translation>
    </message>
    <message>
        <source> degree</source>
        <translation> stupeň</translation>
    </message>
    <message>
        <source> degrees</source>
        <translation>stupne</translation>
    </message>
    <message>
        <source> unknown</source>
        <translation> neznáma</translation>
    </message>
    <message>
        <source>Received %1 of %2 bytes</source>
        <translation type="unfinished">Prijatých %1 z %2 bajtov</translation>
    </message>
    <message>
        <source>Received %1 bytes (total unknown)</source>
        <translation type="unfinished">Prijatých %1 (spolu neznáme)</translation>
    </message>
    <message>
        <source>Not connected</source>
        <translation>Nepripojený</translation>
    </message>
    <message>
        <source>Looking up &apos;%1&apos;</source>
        <translation>Hľadá sa &apos;%1&apos;</translation>
    </message>
    <message>
        <source>Connecting to &apos;%1&apos;</source>
        <translation>Pripája sa k &apos;%1&apos;</translation>
    </message>
    <message>
        <source>Sending request &apos;%1&apos;</source>
        <translation>Posiela sa požiadavka na &apos;%1&apos;</translation>
    </message>
    <message>
        <source>Receiving reply</source>
        <translation>Prijíma sa odpoveď</translation>
    </message>
    <message>
        <source>Response is complete</source>
        <translation>Odpoveď dokončená</translation>
    </message>
    <message>
        <source>Closing down connection</source>
        <translation>Zatvára sa spojenie</translation>
    </message>
    <message>
        <source>Location: </source>
        <comment>

Metadata in GRASS Browser</comment>
        <translation type="obsolete">Lokalita: </translation>
    </message>
    <message>
        <source>&lt;br&gt;Mapset: </source>
        <comment>

Metadata in GRASS Browser</comment>
        <translation type="obsolete">&lt;br&gt;Súbor máp (mapset): </translation>
    </message>
    <message>
        <source>Location: </source>
        <translation>Lokalita: </translation>
    </message>
    <message>
        <source>&lt;br&gt;Mapset: </source>
        <translation>&lt;br&gt;Súbor máp (mapset): </translation>
    </message>
    <message>
        <source>&lt;b&gt;Raster&lt;/b&gt;</source>
        <translation>&lt;b&gt;Raster&lt;/b&gt;</translation>
    </message>
    <message>
        <source>Cannot open raster header</source>
        <translation>Nemožno otvoriť hlavičku rastra</translation>
    </message>
    <message>
        <source>Rows</source>
        <translation>Riadkov</translation>
    </message>
    <message>
        <source>Columns</source>
        <translation>Stĺpcov</translation>
    </message>
    <message>
        <source>N-S resolution</source>
        <translation>Rozlíšenie (S-J)</translation>
    </message>
    <message>
        <source>E-W resolution</source>
        <translation>Rozlíšenie (V-Z)</translation>
    </message>
    <message>
        <source>North</source>
        <translation>Sever</translation>
    </message>
    <message>
        <source>South</source>
        <translation>Juh</translation>
    </message>
    <message>
        <source>East</source>
        <translation>Východ</translation>
    </message>
    <message>
        <source>West</source>
        <translation>Západ</translation>
    </message>
    <message>
        <source>Format</source>
        <translation>Formát</translation>
    </message>
    <message>
        <source>Minimum value</source>
        <translation>Minimálna hodnota</translation>
    </message>
    <message>
        <source>Maximum value</source>
        <translation>Maximálna hodnota</translation>
    </message>
    <message>
        <source>Data source</source>
        <translation>Zdroj údajov</translation>
    </message>
    <message>
        <source>Data description</source>
        <translation>Popis údajov</translation>
    </message>
    <message>
        <source>Comments</source>
        <translation>Komentáre</translation>
    </message>
    <message>
        <source>Categories</source>
        <translation>Kategórie</translation>
    </message>
    <message>
        <source>&lt;b&gt;Vector&lt;/b&gt;</source>
        <translation>&lt;b&gt;Vektor&lt;/b&gt;</translation>
    </message>
    <message>
        <source>Points</source>
        <translation>Body</translation>
    </message>
    <message>
        <source>Lines</source>
        <translation>Línie</translation>
    </message>
    <message>
        <source>Boundaries</source>
        <translation>Hranice</translation>
    </message>
    <message>
        <source>Centroids</source>
        <translation>Centroidy</translation>
    </message>
    <message>
        <source>Faces</source>
        <translation type="unfinished">Faces</translation>
    </message>
    <message>
        <source>Kernels</source>
        <translation type="unfinished">Kernels</translation>
    </message>
    <message>
        <source>Areas</source>
        <translation>Oblasti</translation>
    </message>
    <message>
        <source>Islands</source>
        <translation>Ostrovy</translation>
    </message>
    <message>
        <source>Top</source>
        <translation type="unfinished">Vrch</translation>
    </message>
    <message>
        <source>Bottom</source>
        <translation type="unfinished">Spodok</translation>
    </message>
    <message>
        <source>yes</source>
        <translation>áno</translation>
    </message>
    <message>
        <source>no</source>
        <translation>nie</translation>
    </message>
    <message>
        <source>History&lt;br&gt;</source>
        <translation>História&lt;br&gt;</translation>
    </message>
    <message>
        <source>&lt;b&gt;Layer&lt;/b&gt;</source>
        <translation>&lt;b&gt;Vrstva&lt;/b&gt;</translation>
    </message>
    <message>
        <source>Features</source>
        <translation>Objekty</translation>
    </message>
    <message>
        <source>Driver</source>
        <translation>Ovládač</translation>
    </message>
    <message>
        <source>Database</source>
        <translation>Databáza</translation>
    </message>
    <message>
        <source>Table</source>
        <translation>Tabuľka</translation>
    </message>
    <message>
        <source>Key column</source>
        <translation>Stĺpec s kľúčom</translation>
    </message>
    <message>
        <source>GISBASE is not set.</source>
        <translation>GISBASE nie je nastavená.</translation>
    </message>
    <message>
        <source> is not a GRASS mapset.</source>
        <translation> nie je súbor máp (mapset) GRASSu.</translation>
    </message>
    <message>
        <source>Cannot start </source>
        <translation>Nemožno spustiť </translation>
    </message>
    <message>
        <source>Mapset is already in use.</source>
        <translation>Súbor máp (mapset) sa už používa.</translation>
    </message>
    <message>
        <source>Temporary directory </source>
        <translation>Dočasný adresár </translation>
    </message>
    <message>
        <source> exist but is not writable</source>
        <translation> existuje, ale nemožno do neho zapisovať</translation>
    </message>
    <message>
        <source>Cannot create temporary directory </source>
        <translation>Nemožno vytvoriť dočasný adresár </translation>
    </message>
    <message>
        <source>Cannot create </source>
        <translation>Nemožno vytvoriť </translation>
    </message>
    <message>
        <source>Cannot remove mapset lock: </source>
        <translation>Nemožno odobrať zámok súboru máp (mapsetu): </translation>
    </message>
    <message>
        <source>Warning</source>
        <translation>Upozornenie</translation>
    </message>
    <message>
        <source>Cannot read raster map region</source>
        <translation>Nemožno prečítať región rastrovej mapy</translation>
    </message>
    <message>
        <source>Cannot read vector map region</source>
        <translation>Nemožno prečítať región vektorovej mapy</translation>
    </message>
    <message>
        <source>Cannot read region</source>
        <translation>Nemožno prečítať región</translation>
    </message>
    <message>
        <source>Unable to open </source>
        <translation>Nemožno otvoriť </translation>
    </message>
    <message>
        <source>2.5D shape type not supported</source>
        <translation type="unfinished">Typ 2.5D shape nie je podporovaný</translation>
    </message>
    <message>
        <source>Adding features to 2.5D shapetypes is not supported yet</source>
        <translation type="unfinished">Pridávanie objektov do typu 2.5D shape zatiaľ nie je podporované</translation>
    </message>
    <message>
        <source>Layer cannot be added to</source>
        <translation type="unfinished">Vrstva nemôže byť pridaná</translation>
    </message>
    <message>
        <source>The data provider for this layer does not support the addition of features.</source>
        <translation type="unfinished">Tento nástroj na prístup k údajom nepodporuje pridávanie objektov.</translation>
    </message>
    <message>
        <source>Coordinate transform error</source>
        <translation type="unfinished">Chyba pri tranformácii súradníc</translation>
    </message>
    <message>
        <source>Cannot transform the point to the layers coordinate system</source>
        <translation type="unfinished">Tento bod nemožno transformovať do súradnicového systému vrstiev</translation>
    </message>
    <message>
        <source>Cannot add feature. Unknown WKB type</source>
        <translation type="unfinished">Nemožno pridať objekt. Neznámy typ WKB</translation>
    </message>
    <message>
        <source>Error, could not add island</source>
        <translation type="unfinished">Chyba, nemožno pridať ostrov</translation>
    </message>
    <message>
        <source>A problem with geometry type occured</source>
        <translation type="unfinished">Vyskytol sa problém s typom geometrie</translation>
    </message>
    <message>
        <source>The inserted Ring is not closed</source>
        <translation type="unfinished">Vkladaný prstenec nie je uzavretý</translation>
    </message>
    <message>
        <source>The inserted Ring is not a valid geometry</source>
        <translation type="unfinished">Vkladaný prstenec nemá platnú geometriu</translation>
    </message>
    <message>
        <source>The inserted Ring crosses existing rings</source>
        <translation type="unfinished">Vkladaný prstenec sa križuje s existujúcimi prstencami</translation>
    </message>
    <message>
        <source>The inserted Ring is not contained in a feature</source>
        <translation type="unfinished">Vkladaný prstenec nie je obsiahnutý v objekte</translation>
    </message>
    <message>
        <source>An unknown error occured</source>
        <translation type="unfinished">Vyskytla sa neznáma chyba</translation>
    </message>
    <message>
        <source>Error, could not add ring</source>
        <translation type="unfinished">Chyba, nemožno pridať prstenec</translation>
    </message>
    <message>
        <source>To select features, you must choose a vector layer by clicking on its name in the legend</source>
        <translation type="unfinished">Na to, aby bolo možné vybrať objekty, je potrebné najskôr vybrať aktívnu vrstvu kliknutím na jej meno v okne Legenda</translation>
    </message>
    <message>
        <source>Python error</source>
        <translation type="unfinished">Chyba Pythonu</translation>
    </message>
    <message>
        <source>Couldn&apos;t load plugin </source>
        <translation type="unfinished">Nemožno nahrať zásuvný modul </translation>
    </message>
    <message>
        <source> due an error when calling its classFactory() method</source>
        <translation type="unfinished"> kvôli chybe pri volaní metódy z classFactory()</translation>
    </message>
    <message>
        <source> due an error when calling its initGui() method</source>
        <translation type="unfinished"> kvôli chybe pri volaní metódy initGui()</translation>
    </message>
    <message>
        <source>Error while unloading plugin </source>
        <translation type="unfinished">Chyba pri odhrávaní zásuvného modulu z pamäte </translation>
    </message>
    <message>
        <source>Regular expressions on numeric values don&apos;t make sense. Use comparison instead.</source>
        <translation type="unfinished">Regulárne výrazy nemajú pri numerických hodnotách zmysel. Namiesto toho použite porovnávanie hodnôt.</translation>
    </message>
    <message>
        <source>Geoprocessing functions for working with PostgreSQL/PostGIS layers</source>
        <translation type="unfinished">Funkcie geoprocessing-u sú dostupné pre vrstvy PostgreSQL/PostGIS</translation>
    </message>
    <message>
        <source>Where is &apos;</source>
        <translation type="unfinished">Kde je &apos;</translation>
    </message>
    <message>
        <source>original location: </source>
        <translation type="unfinished">pôvodné umiestnenie: </translation>
    </message>
    <message>
        <source>To identify features, you must choose an active layer by clicking on its name in the legend</source>
        <translation type="unfinished">Na to, aby bolo možné identifikovať objekty, je potrebné najskôr vybrať aktívnu vrstvu kliknutím na jej meno v okne Legenda</translation>
    </message>
    <message>
        <source>PostgreSQL Geoprocessing</source>
        <translation type="unfinished">PostgreSQL Geoprocessing</translation>
    </message>
    <message>
        <source>Quick Print</source>
        <translation>Rýchla tlač</translation>
    </message>
    <message>
        <source>Quick Print is a plugin to quickly print a map with minimal effort.</source>
        <translation type="unfinished">Rýchla tlač je zásuvný modul na rýchle vytlačenie mapy s minimálnym úsilím.</translation>
    </message>
    <message>
        <source>Could not remove polygon intersection</source>
        <translation type="unfinished">Nemožno odobrať prienik polygónu</translation>
    </message>
    <message>
        <source>Currently only filebased datasets are supported</source>
        <translation type="obsolete">V súčasnosti sú podporované len súbory údajov na základe súboru</translation>
    </message>
    <message>
        <source>Loaded default style file from </source>
        <translation type="obsolete">Nahratý predvolený štýl súboru z </translation>
    </message>
    <message>
        <source>The directory containing your dataset needs to be writeable!</source>
        <translation type="unfinished">Adresár obsahujúci vašu súpravu údajov musí byť zapisovateľný!</translation>
    </message>
    <message>
        <source>Created default style file as </source>
        <translation type="unfinished">Vytvorený predvolený štýl súboru ako </translation>
    </message>
    <message>
        <source>ERROR: Failed to created default style file as </source>
        <translation type="obsolete">CHYBA: Vytvorenie súboru s predvoleným štýlom zlyhalo </translation>
    </message>
    <message>
        <source>File could not been opened.</source>
        <translation type="obsolete">Súbor nemožno otvoriť.</translation>
    </message>
    <message>
        <source>Couldn&apos;t load SIP module.</source>
        <translation>Nemožno nahrať modul SIP.</translation>
    </message>
    <message>
        <source>Python support will be disabled.</source>
        <translation type="unfinished">Podpora Pythonu bude vypnutá.</translation>
    </message>
    <message>
        <source>Couldn&apos;t load PyQt4.</source>
        <translation>Nemožno nahrať PyQt4.</translation>
    </message>
    <message>
        <source>Couldn&apos;t load PyQGIS.</source>
        <translation>Nemožno nahrať PyQGIS.</translation>
    </message>
    <message>
        <source>An error has occured while executing Python code:</source>
        <translation>Pri súšťaní kódu v Pythone nastala chyba:</translation>
    </message>
    <message>
        <source>Python version:</source>
        <translation>Verzia Pythonu:</translation>
    </message>
    <message>
        <source>Python path:</source>
        <translation>Cesta k Pythonu:</translation>
    </message>
    <message>
        <source>An error occured during execution of following code:</source>
        <translation>Pri spúšťaní nasledujúceho kódu nastala chyba:</translation>
    </message>
    <message>
        <source> is not writeable.</source>
        <translation>nie je zapisovateľný.</translation>
    </message>
    <message>
        <source>Please adjust permissions (if possible) and try again.</source>
        <translation type="unfinished">Pridajte súboru práva (ak je to možné) a skúste znova.</translation>
    </message>
    <message>
        <source>Uncatched fatal GRASS error</source>
        <translation type="unfinished">Neodchytená fatálna chyba GRASSu</translation>
    </message>
    <message>
        <source>ERROR: Failed to created default style file as %1 Check file permissions and retry.</source>
        <translation type="unfinished">CHYBA: Zlyhalo vytvorenie súboru %1 s predvoleným štýlom Skontrolujte práva súboru a skúste znova.</translation>
    </message>
    <message>
        <source>Coordinate Capture</source>
        <translation>Odchyt súradníc</translation>
    </message>
    <message>
        <source>Capture mouse coordinates in different CRS</source>
        <translation type="unfinished">Odchytí (zaznamená) súradnice myši v rôznych referenčných súradnicových systémoch</translation>
    </message>
    <message>
        <source>Location: </source>
        <comment>Metadata in GRASS Browser</comment>
        <translation type="unfinished">Lokalita: </translation>
    </message>
    <message>
        <source>&lt;br&gt;Mapset: </source>
        <comment>Metadata in GRASS Browser</comment>
        <translation type="unfinished">&lt;br&gt;Súbor máp (mapset): </translation>
    </message>
    <message>
        <source>CRS Exception</source>
        <translation type="unfinished">Výnimka CRS</translation>
    </message>
    <message>
        <source>Selection extends beyond layer&apos;s coordinate system.</source>
        <translation type="unfinished">Výber presahuje mimo súradnicového systému vrstvy.</translation>
    </message>
    <message>
        <source>Legend</source>
        <translation type="unfinished">Legenda</translation>
    </message>
    <message>
        <source>Dxf2Shp Converter</source>
        <translation>Prevodník Dxf2Shp</translation>
    </message>
    <message>
        <source>Converts from dxf to shp file format</source>
        <translation type="unfinished">Prevádza z formátu dxf do shp</translation>
    </message>
    <message>
        <source>Interpolating...</source>
        <translation>Interpoluje sa...</translation>
    </message>
    <message>
        <source>Abort</source>
        <translation>Prerušiť</translation>
    </message>
    <message>
        <source>Interpolation plugin</source>
        <translation>Zásuvný modul na interpoláciu</translation>
    </message>
    <message>
        <source>A plugin for interpolation based on vertices of a vector layer</source>
        <translation>Zásuvný modul slúžiaci na interpoláciu založenú na uzloch vektorovej vrstvy</translation>
    </message>
    <message>
        <source>Version 0.001</source>
        <translation>Verzia 0.001</translation>
    </message>
    <message>
        <source>OGR Layer Converter</source>
        <translation>Prevodník OGR vrstiev</translation>
    </message>
    <message>
        <source>Translates vector layers between formats supported by OGR library</source>
        <translation type="unfinished">Prevádza vektorové vrstvy medzi formátmi podporovanými knižnicou OGR</translation>
    </message>
    <message>
        <source>Loading style file </source>
        <translation type="unfinished">Načítanie súboru štýlu </translation>
    </message>
    <message>
        <source> failed because:</source>
        <translation type="unfinished">zlyhalo kvôli nasledujúcej chybe:</translation>
    </message>
    <message>
        <source>Could not save symbology because:</source>
        <translation type="unfinished">Nemožno uložiť symboliku kvôli:</translation>
    </message>
    <message>
        <source>Unable to save to file. Your project may be corrupted on disk. Try clearing some space on the volume and check file permissions before pressing save again.</source>
        <translation type="unfinished">Nemožno uložiť do súboru. Váš projekt môže byť poškodený na disku. Uvoľnite miesto na disku a skotrolujte práva na zápis pred opätovným kliknutím na tlačidlo Uložiť.</translation>
    </message>
    <message>
        <source>Error Loading Plugin</source>
        <translation type="unfinished">Chyba pri nahrávaní zásuvného modulu</translation>
    </message>
    <message>
        <source>There was an error loading a plugin.The following diagnostic information may help the QGIS developers resolve the issue:
%1.</source>
        <translation type="unfinished">Pri nahrávaní zásuvného modulu sa vyskytla chyba. Nasledujúce diagnostické informácie môžu pomôcť vývojárom QGISu vyriešiť problém:
%1.</translation>
    </message>
    <message>
        <source>Error when reading metadata of plugin </source>
        <translation type="unfinished">Chyba pri čítaní meta údajov zásuvného modulu </translation>
    </message>
</context>
<context>
    <name>QgisApp</name>
    <message>
        <source>Quantum GIS - </source>
        <translation>Quantum GIS -</translation>
    </message>
    <message>
        <source>Version</source>
        <translation>verzia </translation>
    </message>
    <message>
        <source>is not a valid or recognized data source</source>
        <translation>nie je platný alebo rozpoznaný zdroj údajov</translation>
    </message>
    <message>
        <source>Invalid Data Source</source>
        <translation>Chybný zdroj údajov</translation>
    </message>
    <message>
        <source>No Layer Selected</source>
        <translation>Nie je vybratá žiadna vrstva</translation>
    </message>
    <message>
        <source>No MapLayer Plugins</source>
        <translation type="obsolete">Žiadne zásuvné moduly na prácu s vrstvou</translation>
    </message>
    <message>
        <source>No MapLayer plugins in ../plugins/maplayer</source>
        <translation type="obsolete">Žiadne zásuvné moduly na prácu s vrstvou v adresári ../plugins/maplayer</translation>
    </message>
    <message>
        <source>No Plugins</source>
        <translation type="obsolete">Žiadne zásuvné moduly</translation>
    </message>
    <message>
        <source>No plugins found in ../plugins. To test plugins, start qgis from the src directory</source>
        <translation type="obsolete">V adresári ../plugins sa nenašli žiadne zásuvné moduly. Na odskúšanie zásuvných modulov spustite qgis z adresára src</translation>
    </message>
    <message>
        <source>Name</source>
        <translation type="obsolete">Meno</translation>
    </message>
    <message>
        <source>Plugin %1 is named %2</source>
        <translation type="obsolete">Zásuvný modul %1 má názov %2</translation>
    </message>
    <message>
        <source>Plugin Information</source>
        <translation type="obsolete">Informácie o zásuvných moduloch</translation>
    </message>
    <message>
        <source>QGis loaded the following plugin:</source>
        <translation type="obsolete">QGis načítal nasledujúci zásuvný modul:</translation>
    </message>
    <message>
        <source>Name: %1</source>
        <translation type="obsolete">Meno: %1</translation>
    </message>
    <message>
        <source>Version: %1</source>
        <translation type="obsolete">Verzia: %1</translation>
    </message>
    <message>
        <source>Description: %1</source>
        <translation type="obsolete">Popis: %1</translation>
    </message>
    <message>
        <source>Unable to Load Plugin</source>
        <translation type="obsolete">Nie je možné nahrať zásuvný modul</translation>
    </message>
    <message>
        <source>QGIS was unable to load the plugin from: %1</source>
        <translation type="obsolete">QGIS nemohol nahrať zásuvný modul z: %1</translation>
    </message>
    <message>
        <source>There is a new version of QGIS available</source>
        <translation>K dispozícii je nová verzia QGIS</translation>
    </message>
    <message>
        <source>You are running a development version of QGIS</source>
        <translation>Používate vývojovú verziu QGIS</translation>
    </message>
    <message>
        <source>You are running the current version of QGIS</source>
        <translation>Používate aktuálnu verziu QGISu</translation>
    </message>
    <message>
        <source>Would you like more information?</source>
        <translation>Prajete si viac informácií?</translation>
    </message>
    <message>
        <source>QGIS Version Information</source>
        <translation>Informácie o verzii QGIS</translation>
    </message>
    <message>
        <source>Unable to get current version information from server</source>
        <translation>Zo servera nie je možné získať informáciu o aktuálnej verzii</translation>
    </message>
    <message>
        <source>Connection refused - server may be down</source>
        <translation>Spojenie odmietnuté - server je zrejme vypnutý</translation>
    </message>
    <message>
        <source>QGIS server was not found</source>
        <translation>Nebol nájdený server QGIS</translation>
    </message>
    <message>
        <source>Invalid Layer</source>
        <translation>Chybná vrstva</translation>
    </message>
    <message>
        <source>%1 is an invalid layer and cannot be loaded.</source>
        <translation>Vrstva %1 je chybná a nemôže byť nahratá.</translation>
    </message>
    <message>
        <source>Error Loading Plugin</source>
        <translation type="obsolete">Chyba pri nahrávaní zásuvného modulu</translation>
    </message>
    <message>
        <source>There was an error loading %1.</source>
        <translation type="obsolete">Pri nahrávaní %1 sa vyskytla chyba.</translation>
    </message>
    <message>
        <source>Saved map image to</source>
        <translation>Obrázok mapy uložený do</translation>
    </message>
    <message>
        <source>Choose a filename to save the map image as</source>
        <translation type="obsolete">Vyberte meno súboru, do ktorého sa má uložiť obrázok mapy</translation>
    </message>
    <message>
        <source>Extents: </source>
        <translation>Rozsah: </translation>
    </message>
    <message>
        <source>Problem deleting features</source>
        <translation type="unfinished">Problém pri mazaní objektov</translation>
    </message>
    <message>
        <source>A problem occured during deletion of features</source>
        <translation type="unfinished">Pri mazaní objektov sa vyskytol problém</translation>
    </message>
    <message>
        <source>No Vector Layer Selected</source>
        <translation>Nie je vybratá žiadna vrstva</translation>
    </message>
    <message>
        <source>Deleting features only works on vector layers</source>
        <translation type="unfinished">Mazanie objektov funguje iba pri vektorových vrstvách</translation>
    </message>
    <message>
        <source>To delete features, you must select a vector layer in the legend</source>
        <translation type="unfinished">Na vymazanie objektu je treba vybrať vektorovú vrstvu v okne Legenda</translation>
    </message>
    <message>
        <source>Render</source>
        <translation>Vykresľovanie</translation>
    </message>
    <message>
        <source>Choose a QGIS project file</source>
        <translation type="unfinished">Vyberte súbor QGIS projektu</translation>
    </message>
    <message>
        <source>Unable to save project</source>
        <translation>Nie je možné uložiť projekt</translation>
    </message>
    <message>
        <source>Unable to save project to </source>
        <translation type="unfinished">Nie je možné uložiť projekt do </translation>
    </message>
    <message>
        <source>Map legend that displays all the layers currently on the map canvas. Click on the check box to turn a layer on or off. Double click on a layer in the legend to customize its appearance and set other properties.</source>
        <translation>Okno Legenda, ktoré zobrazuje všetky vrstvy nachádzajúce sa na mapovom plátne. Kliknutím na zaškrtávacie políčko sa zapne alebo vypne vrstva. Dvojklikom na vrstvu v okne Legenda možno prispôsobiť jej vzhľad a nastaviť ostatné vlastnosti.</translation>
    </message>
    <message>
        <source>Map overview canvas. This canvas can be used to display a locator map that shows the current extent of the map canvas. The current extent is shown as a red rectangle. Any layer on the map can be added to the overview canvas.</source>
        <translation>Okno Mapový prehľad. Táto oblasť môže byť použitá na zobrazenie polohy mapy, kde vidno práve vyobrazený výrez z mapy. Aktuálne zobrazený rozsah je znázornený červeným obdĺžnikom. Do mapového prehľadu je možné pridať ktorúkoľvek vrstvu.</translation>
    </message>
    <message>
        <source>&amp;Plugins</source>
        <translation>Zásuvné &amp;moduly</translation>
    </message>
    <message>
        <source>Displays the current map scale</source>
        <translation>Zobrazuje aktuálnu mierku mapy</translation>
    </message>
    <message>
        <source>When checked, the map layers are rendered in response to map navigation commands and other events. When not checked, no rendering is done. This allows you to add a large number of layers and symbolize them before rendering.</source>
        <translation>Pokiaľ je políčko zaškrtnuté, mapové vrstvy sú vykresľované v náväznosti na príkazy mapovej navigácie a ďalšie udalosti. Ak nie je zaškrtnuté, nič sa nevykresľuje. To dovoľuje pridať veľké množstvo vrstiev a upraviť ich symboliku ešte pred ich vykreslením.</translation>
    </message>
    <message>
        <source>Open an OGR Supported Vector Layer</source>
        <translation type="unfinished">Otvoriť vektorovú vrstvu podporovanú knižnicou OGR</translation>
    </message>
    <message>
        <source>QGIS Project Read Error</source>
        <translation>Chyba pri čítaní QGIS projektu</translation>
    </message>
    <message>
        <source>Open a GDAL Supported Raster Data Source</source>
        <translation type="unfinished">Otvoriť rastrový zdroj údajov podporovaný knižnicou GDAL</translation>
    </message>
    <message>
        <source>Toggle map rendering</source>
        <translation>Prepínač vykresľovania mapy</translation>
    </message>
    <message>
        <source>This icon shows whether on the fly projection is enabled or not. Click the icon to bring up the project properties dialog to alter this behaviour.</source>
        <translation type="obsolete">Táto ikona ukazuje, či je povolený (zapnutý) priamy prevod medzi mapovými zobrazeniami. Kliknutím na túto ikonu sa otvorí dialóg Vlastnosti projektu, kde je možné zmeniť toto nastavenie.</translation>
    </message>
    <message>
        <source>Projection status - Click to open projection dialog</source>
        <translation type="obsolete">Stav zobrazenia - kliknutím sa otvorí menu mapových zobrazení</translation>
    </message>
    <message>
        <source>Try to find missing layers?</source>
        <translation type="unfinished">Pokúsiť sa nájsť chýbajúce vrstvy?</translation>
    </message>
    <message>
        <source>Save As</source>
        <translation>Uložiť ako</translation>
    </message>
    <message>
        <source>Choose a QGIS project file to open</source>
        <translation type="unfinished">Vyberte súbor QGIS projektu, ktorý chcete otvoriť</translation>
    </message>
    <message>
        <source>Saved project to:</source>
        <translation>Projekt uložený do:</translation>
    </message>
    <message>
        <source>New features</source>
        <translation>Nové vlastnosti</translation>
    </message>
    <message>
        <source>Unable to open project</source>
        <translation>Nie je možné otvoriť projekt</translation>
    </message>
    <message>
        <source>Unable to save project </source>
        <translation>Nie je možné uložiť projekt </translation>
    </message>
    <message>
        <source>QGIS: Unable to load project</source>
        <translation type="unfinished">QGIS: Nie je možné načítať projekt</translation>
    </message>
    <message>
        <source>Unable to load project </source>
        <translation>Nemožno načítať projekt </translation>
    </message>
    <message>
        <source>Layer is not valid</source>
        <translation>Neplatná vrstva</translation>
    </message>
    <message>
        <source>The layer is not a valid layer and can not be added to the map</source>
        <translation type="unfinished">Táto vrstva nie je platnou vrstvou a nemôže byť pridaná do mapy</translation>
    </message>
    <message>
        <source>Save?</source>
        <translation>Uložiť?</translation>
    </message>
    <message>
        <source>Do you want to save the current project?</source>
        <translation>Uložiť aktuálny projekt?</translation>
    </message>
    <message>
        <source>Show all layers</source>
        <translation>Ukázať všetky vrstvy</translation>
    </message>
    <message>
        <source>Hide all layers</source>
        <translation>Skryť všetky vrstvy</translation>
    </message>
    <message>
        <source>Clipboard contents set to: </source>
        <translation type="obsolete">Obsah schránky uložený do: </translation>
    </message>
    <message>
        <source> is not a valid or recognized raster data source</source>
        <translation type="unfinished">nie je platný, alebo rozpoznaný zdroj údajov</translation>
    </message>
    <message>
        <source> is not a supported raster data source</source>
        <translation> nie je podporovaný rastrový zdroj údajov</translation>
    </message>
    <message>
        <source>Unsupported Data Source</source>
        <translation>Nepodporovaný zdroj údajov</translation>
    </message>
    <message>
        <source>New Bookmark</source>
        <translation>Nová záložka</translation>
    </message>
    <message>
        <source>Enter a name for the new bookmark:</source>
        <translation>Zadajte meno novej záložky:</translation>
    </message>
    <message>
        <source>Error</source>
        <translation>Chyba</translation>
    </message>
    <message>
        <source>Unable to create the bookmark. Your user database may be missing or corrupted</source>
        <translation>Nie je možné vytvoriť záložku. Databáza používateľa chýba, alebo je poškodená</translation>
    </message>
    <message>
        <source>New Project</source>
        <translation>Nový projekt</translation>
    </message>
    <message>
        <source>&amp;Open Project...</source>
        <translation>&amp;Otvoriť projekt...</translation>
    </message>
    <message>
        <source>Open a Project</source>
        <translation>Otvoriť projekt</translation>
    </message>
    <message>
        <source>&amp;Save Project</source>
        <translation>&amp;Uložiť projekt</translation>
    </message>
    <message>
        <source>Save Project &amp;As...</source>
        <translation>Uložiť projekt &amp;ako...</translation>
    </message>
    <message>
        <source>Save Project under a new name</source>
        <translation>Uložiť projekt pod novým menom</translation>
    </message>
    <message>
        <source>Print</source>
        <translation type="obsolete">Tlačiť</translation>
    </message>
    <message>
        <source>Save map as image</source>
        <translation>Uložiť mapu ako obrázok</translation>
    </message>
    <message>
        <source>Exit</source>
        <translation>Koniec</translation>
    </message>
    <message>
        <source>Exit QGIS</source>
        <translation>Ukončiť QGIS</translation>
    </message>
    <message>
        <source>Add a Vector Layer</source>
        <translation>Pridať vektorovú vrstvu</translation>
    </message>
    <message>
        <source>Add a Raster Layer</source>
        <translation>Pridať rastrovú vrstvu</translation>
    </message>
    <message>
        <source>Add a PostGIS Layer</source>
        <translation>Pridať vrstvu PostGIS</translation>
    </message>
    <message>
        <source>Remove Layer</source>
        <translation>Odobrať vrstvu</translation>
    </message>
    <message>
        <source>Remove a Layer</source>
        <translation>Odoberie vrstvu</translation>
    </message>
    <message>
        <source>Create a New Vector Layer</source>
        <translation>Vytvoriť novú vektorovú vrstvu</translation>
    </message>
    <message>
        <source>Add All To Overview</source>
        <translation type="obsolete">Pridať všetky do Prehľadu</translation>
    </message>
    <message>
        <source>Show all layers in the overview map</source>
        <translation>Ukázať všetky vrstvy v Prehľade mapy</translation>
    </message>
    <message>
        <source>Remove All From Overview</source>
        <translation>Odobrať všetky z Prehľadu</translation>
    </message>
    <message>
        <source>Remove all layers from overview map</source>
        <translation>Odoberať všetky vrstvy z Prehľadu mapy</translation>
    </message>
    <message>
        <source>Show All Layers</source>
        <translation>Ukázať všetky vrstvy</translation>
    </message>
    <message>
        <source>Hide All Layers</source>
        <translation>Skryť všetky vrstvy</translation>
    </message>
    <message>
        <source>Set project properties</source>
        <translation>Nastaviť vlastnosti projektu</translation>
    </message>
    <message>
        <source>Change various QGIS options</source>
        <translation>Zmeniť rôzne vlastnosti QGIS</translation>
    </message>
    <message>
        <source>Manage custom projections</source>
        <translation type="obsolete">Správa vlastných kartografických zobrazení</translation>
    </message>
    <message>
        <source>Help Contents</source>
        <translation>Obsah Pomocníka</translation>
    </message>
    <message>
        <source>Help Documentation</source>
        <translation>Sprievodná dokumentácia</translation>
    </message>
    <message>
        <source>Qgis Home Page</source>
        <translation type="obsolete">Domovská stránka Qgis</translation>
    </message>
    <message>
        <source>QGIS Home Page</source>
        <translation>Domovská stránka QGIS</translation>
    </message>
    <message>
        <source>About</source>
        <translation>O programe</translation>
    </message>
    <message>
        <source>About QGIS</source>
        <translation>O programe QGIS</translation>
    </message>
    <message>
        <source>Check Qgis Version</source>
        <translation>Skontrolovať verziu Qgis</translation>
    </message>
    <message>
        <source>Check if your QGIS version is up to date (requires internet access)</source>
        <translation>Skontroluje, či je vaša verzia QGIS aktuálna (to vyžaduje prístup k internetu)</translation>
    </message>
    <message>
        <source>Refresh</source>
        <translation>Obnoviť</translation>
    </message>
    <message>
        <source>Refresh Map</source>
        <translation>Obnoviť mapu</translation>
    </message>
    <message>
        <source>Zoom In</source>
        <translation>Priblížiť</translation>
    </message>
    <message>
        <source>Zoom Out</source>
        <translation>Oddialiť</translation>
    </message>
    <message>
        <source>Zoom Full</source>
        <translation>Celá mapa</translation>
    </message>
    <message>
        <source>Zoom to Full Extents</source>
        <translation>Zmeniť pohľad na veľkosť celej mapy</translation>
    </message>
    <message>
        <source>Zoom To Selection</source>
        <translation type="obsolete">Na veľkosť výberu</translation>
    </message>
    <message>
        <source>Zoom to selection</source>
        <translation type="obsolete">Zmení pohľad na veľkosť výberu</translation>
    </message>
    <message>
        <source>Pan Map</source>
        <translation>Posun mapy</translation>
    </message>
    <message>
        <source>Pan the map</source>
        <translation>Posunúť mapu</translation>
    </message>
    <message>
        <source>Zoom Last</source>
        <translation>Predchádzajúci pohľad</translation>
    </message>
    <message>
        <source>Zoom to Last Extent</source>
        <translation>Predchádzajúci pohľad</translation>
    </message>
    <message>
        <source>Zoom To Layer</source>
        <translation type="obsolete">Na veľkosť vrstvy</translation>
    </message>
    <message>
        <source>Zoom to Layer</source>
        <translation>Zmení pohľad na veľkosť vrstvy</translation>
    </message>
    <message>
        <source>Identify Features</source>
        <translation>Identifikovať objekty</translation>
    </message>
    <message>
        <source>Click on features to identify them</source>
        <translation>Kliknutím na objekty ich možno identifikovať</translation>
    </message>
    <message>
        <source>Select Features</source>
        <translation>Vybrať objekty</translation>
    </message>
    <message>
        <source>Open Table</source>
        <translation type="obsolete">Otvoriť tabuľku</translation>
    </message>
    <message>
        <source>Measure Line </source>
        <translation> Merať vzdialenosť</translation>
    </message>
    <message>
        <source>Measure a Line</source>
        <translation>Meranie vzdialeností</translation>
    </message>
    <message>
        <source>Measure Area</source>
        <translation>Meranie rozlohy</translation>
    </message>
    <message>
        <source>Measure an Area</source>
        <translation>Merať rozlohu</translation>
    </message>
    <message>
        <source>Show Bookmarks</source>
        <translation>Ukázať záložky</translation>
    </message>
    <message>
        <source>Add Web Mapping Server Layer</source>
        <translation type="obsolete">Pridať vrstvu z Web Mapping Server-a</translation>
    </message>
    <message>
        <source>In Overview</source>
        <translation type="obsolete">Do Prehľadu</translation>
    </message>
    <message>
        <source>Add current layer to overview map</source>
        <translation>Pridať aktívnu vrstvu do Prehľadu mapy</translation>
    </message>
    <message>
        <source>Open the plugin manager</source>
        <translation>Otvoriť správcu zásuvných modulov</translation>
    </message>
    <message>
        <source>Capture Point</source>
        <translation>Získať bod</translation>
    </message>
    <message>
        <source>Capture Points</source>
        <translation>Získať body</translation>
    </message>
    <message>
        <source>Capture Line</source>
        <translation>Získať líniu</translation>
    </message>
    <message>
        <source>Capture Lines</source>
        <translation>Získať línie</translation>
    </message>
    <message>
        <source>Capture Polygon</source>
        <translation>Získať polygón</translation>
    </message>
    <message>
        <source>Capture Polygons</source>
        <translation>Získať polygóny</translation>
    </message>
    <message>
        <source>&amp;File</source>
        <translation>&amp;Súbor</translation>
    </message>
    <message>
        <source>&amp;View</source>
        <translation>Pohľa&amp;d</translation>
    </message>
    <message>
        <source>&amp;Layer</source>
        <translation>&amp;Vrstva</translation>
    </message>
    <message>
        <source>&amp;Settings</source>
        <translation>&amp;Nastavenia</translation>
    </message>
    <message>
        <source>&amp;Help</source>
        <translation>&amp;Pomocník</translation>
    </message>
    <message>
        <source>File</source>
        <translation>Súbor</translation>
    </message>
    <message>
        <source>Manage Layers</source>
        <translation>Správa vrstiev</translation>
    </message>
    <message>
        <source>Help</source>
        <translation>Pomocník</translation>
    </message>
    <message>
        <source>Digitizing</source>
        <translation>Digitalizácia</translation>
    </message>
    <message>
        <source>Map Navigation</source>
        <translation>Navigácia na mape</translation>
    </message>
    <message>
        <source>Attributes</source>
        <translation>Atribúty</translation>
    </message>
    <message>
        <source>Plugins</source>
        <translation>Zásuvné moduly</translation>
    </message>
    <message>
        <source>Ready</source>
        <translation>Pripravený</translation>
    </message>
    <message>
        <source>Delete Selected</source>
        <translation>Zmazať vybrané</translation>
    </message>
    <message>
        <source>Add Vertex</source>
        <translation>Pridať uzol</translation>
    </message>
    <message>
        <source>Delete Vertex</source>
        <translation>Vymazať uzol</translation>
    </message>
    <message>
        <source>Move Vertex</source>
        <translation>Premiestniť uzol</translation>
    </message>
    <message>
        <source>QGIS - Changes in SVN Since Last Release</source>
        <translation>QGIS - Zmeny v SVN od ostatného vydania</translation>
    </message>
    <message>
        <source>&amp;New Project</source>
        <translation>&amp;Nový projekt</translation>
    </message>
    <message>
        <source>&amp;Print...</source>
        <translation type="obsolete">&amp;Tlačiť...</translation>
    </message>
    <message>
        <source>Save as Image...</source>
        <translation>Uložiť ako obrázok...</translation>
    </message>
    <message>
        <source>Add a Vector Layer...</source>
        <translation type="obsolete">Pridať vektorovú vrstvu...</translation>
    </message>
    <message>
        <source>Add a Raster Layer...</source>
        <translation type="obsolete">Pridať rastrovú vrstvu...</translation>
    </message>
    <message>
        <source>Add a PostGIS Layer...</source>
        <translation type="obsolete">Pridať vrstvu PostGIS...</translation>
    </message>
    <message>
        <source>New Vector Layer...</source>
        <translation>Nová vektorová vrstva...</translation>
    </message>
    <message>
        <source>Project Properties...</source>
        <translation>Vlastnosti projektu...</translation>
    </message>
    <message>
        <source>Options...</source>
        <translation>Vlastnosti...</translation>
    </message>
    <message>
        <source>Custom Projection...</source>
        <translation type="obsolete">Vlastné zobrazenia...</translation>
    </message>
    <message>
        <source>New Bookmark...</source>
        <translation>Nová záložka...</translation>
    </message>
    <message>
        <source>Add WMS Layer...</source>
        <translation>Pridať vrstvu WMS...</translation>
    </message>
    <message>
        <source>Plugin Manager...</source>
        <translation type="obsolete">Správca zásuvných modulov...</translation>
    </message>
    <message>
        <source>Reading settings</source>
        <translation>Čítajú sa nastavenia</translation>
    </message>
    <message>
        <source>Setting up the GUI</source>
        <translation>Nastavuje sa grafické používateľské rozhranie</translation>
    </message>
    <message>
        <source>Checking database</source>
        <translation>Kontroluje sa databáza</translation>
    </message>
    <message>
        <source>Restoring loaded plugins</source>
        <translation>Obnovujú sa nahraté zásuvné moduly</translation>
    </message>
    <message>
        <source>Initializing file filters</source>
        <translation>Inicializujú sa filtre súborov</translation>
    </message>
    <message>
        <source>Restoring window state</source>
        <translation>Obnovuje sa stav okna</translation>
    </message>
    <message>
        <source>QGIS Ready!</source>
        <translation>QGIS pripravený!</translation>
    </message>
    <message>
        <source>&amp;Open Recent Projects</source>
        <translation>Nedávno &amp;otvorené projekty</translation>
    </message>
    <message>
        <source>Choose a filename to save the QGIS project file as</source>
        <translation type="obsolete">Vyberte meno súboru do ktorého sa uloží súbor projektu QGIS</translation>
    </message>
    <message>
        <source>Ctrl+N</source>
        <comment>New Project</comment>
        <translation type="obsolete">Ctrl+N</translation>
    </message>
    <message>
        <source>Ctrl+O</source>
        <comment>Open a Project</comment>
        <translation type="obsolete">Ctrl+O</translation>
    </message>
    <message>
        <source>Ctrl+S</source>
        <comment>Save Project</comment>
        <translation type="obsolete">Ctrl+S</translation>
    </message>
    <message>
        <source>Save Project</source>
        <translation>Uložiť projekt</translation>
    </message>
    <message>
        <source>Ctrl+A</source>
        <comment>Save Project under a new name</comment>
        <translation type="obsolete">Ctrl+A</translation>
    </message>
    <message>
        <source>Ctrl+P</source>
        <comment>Print</comment>
        <translation type="obsolete">Ctrl+P</translation>
    </message>
    <message>
        <source>Ctrl+I</source>
        <comment>Save map as image</comment>
        <translation type="obsolete">Ctrl+I</translation>
    </message>
    <message>
        <source>Ctrl+Q</source>
        <comment>Exit QGIS</comment>
        <translation type="obsolete">Ctrl+Q</translation>
    </message>
    <message>
        <source>V</source>
        <comment>Add a Vector Layer</comment>
        <translation type="obsolete">V</translation>
    </message>
    <message>
        <source>R</source>
        <comment>Add a Raster Layer</comment>
        <translation type="obsolete">R</translation>
    </message>
    <message>
        <source>D</source>
        <comment>Add a PostGIS Layer</comment>
        <translation type="obsolete">D</translation>
    </message>
    <message>
        <source>N</source>
        <comment>Create a New Vector Layer</comment>
        <translation type="obsolete">N</translation>
    </message>
    <message>
        <source>Ctrl+D</source>
        <comment>Remove a Layer</comment>
        <translation type="obsolete">Ctrl+D</translation>
    </message>
    <message>
        <source>+</source>
        <comment>Show all layers in the overview map</comment>
        <translation type="obsolete">+</translation>
    </message>
    <message>
        <source>-</source>
        <comment>Remove all layers from overview map</comment>
        <translation type="obsolete">-</translation>
    </message>
    <message>
        <source>S</source>
        <comment>Show all layers</comment>
        <translation type="obsolete">S</translation>
    </message>
    <message>
        <source>H</source>
        <comment>Hide all layers</comment>
        <translation type="obsolete">H</translation>
    </message>
    <message>
        <source>P</source>
        <comment>Set project properties</comment>
        <translation type="obsolete">P</translation>
    </message>
    <message>
        <source>F1</source>
        <comment>Help Documentation</comment>
        <translation type="obsolete">F1</translation>
    </message>
    <message>
        <source>Ctrl+H</source>
        <comment>QGIS Home Page</comment>
        <translation type="obsolete">Ctrl+H</translation>
    </message>
    <message>
        <source>Ctrl+R</source>
        <comment>Refresh Map</comment>
        <translation type="obsolete">Ctrl+R</translation>
    </message>
    <message>
        <source>Ctrl++</source>
        <comment>Zoom In</comment>
        <translation type="obsolete">Ctrl++</translation>
    </message>
    <message>
        <source>Ctrl+-</source>
        <comment>Zoom Out</comment>
        <translation type="obsolete">Ctrl+-</translation>
    </message>
    <message>
        <source>F</source>
        <comment>Zoom to Full Extents</comment>
        <translation type="obsolete">F</translation>
    </message>
    <message>
        <source>Ctrl+F</source>
        <comment>Zoom to selection</comment>
        <translation type="obsolete">Ctrl+F</translation>
    </message>
    <message>
        <source>I</source>
        <comment>Click on features to identify them</comment>
        <translation type="obsolete">I</translation>
    </message>
    <message>
        <source>Ctrl+M</source>
        <comment>Measure a Line</comment>
        <translation type="obsolete">Ctrl+M</translation>
    </message>
    <message>
        <source>Ctrl+J</source>
        <comment>Measure an Area</comment>
        <translation type="obsolete">Ctrl+J</translation>
    </message>
    <message>
        <source>B</source>
        <comment>Show Bookmarks</comment>
        <translation type="obsolete">B</translation>
    </message>
    <message>
        <source>Ctrl+B</source>
        <comment>New Bookmark</comment>
        <translation type="obsolete">Ctrl+B</translation>
    </message>
    <message>
        <source>W</source>
        <comment>Add Web Mapping Server Layer</comment>
        <translation type="obsolete">W</translation>
    </message>
    <message>
        <source>O</source>
        <comment>Add current layer to overview map</comment>
        <translation type="obsolete">O</translation>
    </message>
    <message>
        <source>.</source>
        <comment>Capture Points</comment>
        <translation type="obsolete">.</translation>
    </message>
    <message>
        <source>/</source>
        <comment>Capture Lines</comment>
        <translation type="obsolete">/</translation>
    </message>
    <message>
        <source>Ctrl+/</source>
        <comment>Capture Polygons</comment>
        <translation type="obsolete">Ctrl+/</translation>
    </message>
    <message>
        <source>Ctrl+?</source>
        <comment>Help Documentation (Mac)</comment>
        <translation type="obsolete">Ctrl+?</translation>
    </message>
    <message>
        <source>Cut Features</source>
        <translation>Vystrihnúť objekty</translation>
    </message>
    <message>
        <source>Cut selected features</source>
        <translation>Vystrihnúť vybrané objekty</translation>
    </message>
    <message>
        <source>Copy Features</source>
        <translation>Kopírovať objekty</translation>
    </message>
    <message>
        <source>Copy selected features</source>
        <translation>Kopírovať vybrané objekty</translation>
    </message>
    <message>
        <source>Paste Features</source>
        <translation>Vložiť objekty</translation>
    </message>
    <message>
        <source>Paste selected features</source>
        <translation>Vloží vybrané objekty</translation>
    </message>
    <message>
        <source>Show most toolbars</source>
        <translation type="obsolete">Zobraziť väčšinu panelov s nástrojmi</translation>
    </message>
    <message>
        <source>Hide most toolbars</source>
        <translation type="obsolete">Skryť väčšinu panelov s nástrojmi</translation>
    </message>
    <message>
        <source>Network error while communicating with server</source>
        <translation type="unfinished">Chyba siete pri komunikácii so serverom</translation>
    </message>
    <message>
        <source>Unknown network socket error</source>
        <translation>Neznáma chyba sieťového spojenia</translation>
    </message>
    <message>
        <source>Unable to communicate with QGIS Version server</source>
        <translation>Nemožno nadviazať spojenie so serverom &apos;QGIS Version server&apos;</translation>
    </message>
    <message>
        <source>T</source>
        <comment>Show most toolbars</comment>
        <translation type="obsolete">T</translation>
    </message>
    <message>
        <source>Ctrl+T</source>
        <comment>Hide most toolbars</comment>
        <translation type="obsolete">Ctrl+T</translation>
    </message>
    <message>
        <source>Current map scale</source>
        <translation>Aktuálna mierka mapy</translation>
    </message>
    <message>
        <source>Map coordinates at mouse cursor position</source>
        <translation>Mapové súradnice polohy kurzora myši</translation>
    </message>
    <message>
        <source>Checking provider plugins</source>
        <translation>Kontrolujú sa zásuvné moduly na prístup k údajom</translation>
    </message>
    <message>
        <source>Starting Python</source>
        <translation>Spúšťa sa Python</translation>
    </message>
    <message>
        <source>Toggle editing</source>
        <translation>Prepnúť na úpravy</translation>
    </message>
    <message>
        <source>Toggles the editing state of the current layer</source>
        <translation type="unfinished">Prepne na editáciu aktuálnej vrstvy</translation>
    </message>
    <message>
        <source>Add Ring</source>
        <translation>Pridať prstenec</translation>
    </message>
    <message>
        <source>Add Island</source>
        <translation>Pridať ostrov</translation>
    </message>
    <message>
        <source>Add Island to multipolygon</source>
        <translation>Pridať ostrov do multipolygónu</translation>
    </message>
    <message>
        <source>Python console</source>
        <translation type="obsolete">Konzola Pythonu</translation>
    </message>
    <message>
        <source>Toolbar Visibility...</source>
        <translation type="obsolete">Viditeľnosť panelov...</translation>
    </message>
    <message>
        <source>Scale </source>
        <translation>Mierka </translation>
    </message>
    <message>
        <source>Current map scale (formatted as x:y)</source>
        <translation>Aktuálna mierka mapy (v tvare x:y)</translation>
    </message>
    <message>
        <source>Python error</source>
        <translation type="obsolete">Chyba Pythonu</translation>
    </message>
    <message>
        <source>Error when reading metadata of plugin </source>
        <translation type="obsolete">Chyba pri čítaní meta údajov zásuvného modulu </translation>
    </message>
    <message>
        <source>Provider does not support deletion</source>
        <translation type="unfinished">Správca údajov nepodporuje mazanie</translation>
    </message>
    <message>
        <source>Data provider does not support deleting features</source>
        <translation type="unfinished">Správca údajov nepodporuje mazanie objektov</translation>
    </message>
    <message>
        <source>Layer not editable</source>
        <translation>Vrstva nie je upravovateľná</translation>
    </message>
    <message>
        <source>The current layer is not editable. Choose &apos;Start editing&apos; in the digitizing toolbar.</source>
        <translation type="unfinished">Aktuálna vrstva nie je upravovateľná. Na paneli Digitalizácia kliknite na ikonu &apos;Začať úpravy&apos;.</translation>
    </message>
    <message>
        <source>Invalid scale</source>
        <translation>Neplatná mierka</translation>
    </message>
    <message>
        <source>Move Feature</source>
        <translation>Posun objektu</translation>
    </message>
    <message>
        <source>Split Features</source>
        <translation>Rozdeliť objekty</translation>
    </message>
    <message>
        <source>Map Tips</source>
        <translation>Mapové tipy</translation>
    </message>
    <message>
        <source>Show information about a feature when the mouse is hovered over it</source>
        <translation type="unfinished">Ukáže informáciu o objekte pri prechádzaní myšou ponad objekt</translation>
    </message>
    <message>
        <source>Project file is older</source>
        <translation>Súbor projektu je starší</translation>
    </message>
    <message>
        <source>&lt;p&gt;This project file was saved by an older version of QGIS.</source>
        <translation>&lt;p&gt;Tento súbor projektu bol uložený staršou verziou QGISu.</translation>
    </message>
    <message>
        <source> When saving this project file, QGIS will update it to the latest version, possibly rendering it useless for older versions of QGIS.</source>
        <translation type="unfinished"> Akonáhle uložíte tento súbor projektu, QGIS ho aktualizuje podľa najnovšej verzie, čo môže znamenať, že v starších verziách prestane fungovať.</translation>
    </message>
    <message>
        <source>&lt;p&gt;Even though QGIS developers try to maintain backwards compatibility, some of the information from the old project file might be lost.</source>
        <translation type="unfinished">&lt;p&gt;Aj keď sa vývojarí QGISu snažia zachovať spätnú kompatibilitu, niektoré informácie zo starších súborov projektu sa môžu stratiť.</translation>
    </message>
    <message>
        <source> To improve the quality of QGIS, we appreciate if you file a bug report at %3.</source>
        <translation type="unfinished"> Aby sa zlepšila kvalita QGISu, oceníme ak zaznamenáte správu o chybe na stránke %3</translation>
    </message>
    <message>
        <source> Be sure to include the old project file, and state the version of QGIS you used to discover the error.</source>
        <translation type="unfinished"> Nezabudnite pripojiť starý súbor projektu a stave verzie QGISu, v ktorej sa objeavila chyba.</translation>
    </message>
    <message>
        <source>&lt;p&gt;To remove this warning when opening an older project file, uncheck the box &apos;%5&apos; in the %4 menu.</source>
        <translation type="unfinished">&lt;p&gt;Aby sa nezobrazovalo toto hlásenie v budúcnosti pri otvorení staršieho súboru s projektom, odškrtnite &apos;%5&apos; v %4 menu.</translation>
    </message>
    <message>
        <source>&lt;p&gt;Version of the project file: %1&lt;br&gt;Current version of QGIS: %2</source>
        <translation type="unfinished">&lt;p&gt;Verzia súboru projektu: %1&lt;br&gt;Aktuálna verzia QGISu: %2</translation>
    </message>
    <message>
        <source>&lt;tt&gt;Settings:Options:General&lt;/tt&gt;</source>
        <comment>

Menu path to setting options</comment>
        <translation type="obsolete">&lt;tt&gt;Nastavenia:Vlastnosti:Všeobecné&lt;/tt&gt;</translation>
    </message>
    <message>
        <source>Warn me when opening a project file saved with an older version of QGIS</source>
        <translation type="unfinished">Upozorniť pri otváraní súboru projektu uloženého staršou verziou QGISu</translation>
    </message>
    <message>
        <source>Toggle full screen mode</source>
        <translation type="obsolete">Prepnúč do celoobrazovkového režimu</translation>
    </message>
    <message>
        <source>Ctrl-F</source>
        <comment>Toggle fullscreen mode</comment>
        <translation type="obsolete">Ctrl-F</translation>
    </message>
    <message>
        <source>Toggle fullscreen mode</source>
        <translation>Zapnúť celoobrazovkový režim</translation>
    </message>
    <message>
        <source>Imrovements to digitising capabilities.</source>
        <translation type="obsolete">Vylepšenie digitalizácie.</translation>
    </message>
    <message>
        <source>Supporting default and defined styles (.qml) files for file based vector layers. With styles you can save the symbolisation and other settings associated with a vector layer and they will be loaded whenever you load that layer.</source>
        <translation type="obsolete">Podporuje predvolené a definované súbormi so štýlmi (.qml) pre vektorové vrstvy. Štýli možno použiť na uloženie symboliky a ďalších nastavení prislúchajúcich vektorovej vrstve a možno ho načítať kedykoľvek budete pracovať s danou vektorovou vrstvou.</translation>
    </message>
    <message>
        <source>Improved support for transparency and contrast stretching in raster layers. Support for color ramps in raster layers. Support for non-north up rasters. Many other raster improvements &apos;under the hood&apos;.</source>
        <translation type="obsolete">Vylepšený podpora priehľadnosti a úprava kontrastu v rastrových vrstvách. V rastrových vrstvách podpora pre color ramps. Podpora pre rastre neorientované na sever. Mnoho iných zlepšení &apos;pod pokrievkou&apos;.</translation>
    </message>
    <message>
        <source>Progress bar that displays the statusof rendering layers and other time-intensive operations</source>
        <translation type="obsolete">Ukazovateľ priebehu zobrazuje stav vykresľovania vrstiev a iné oparácie citlivé na čas</translation>
    </message>
    <message>
        <source>Shows the map coordinates at thecurrent cursor position. The display is continuously updated as the mouse is moved.</source>
        <translation type="obsolete">Zobrazuje mapové súradnice aktuálnej polohy kurzora. Zobrazenie je pri pohybe myšiu priebežne aktualizované.</translation>
    </message>
    <message>
        <source>Resource Location Error</source>
        <translation type="unfinished">Chyba pri lokalizovaní zdrojov</translation>
    </message>
    <message>
        <source>Error reading icon resources from: 
 %1
 Quitting...</source>
        <translation type="unfinished">Chyba pri načítavaní zdrojov ikôn z: 
  %1
 Ukončuje sa...</translation>
    </message>
    <message>
        <source>Map canvas. This is where raster and vectorlayers are displayed when added to the map</source>
        <translation type="obsolete">Mapové plátno - je to miesto kde sa zobrazujú rastrové a vektorové vrstvy pridané do mapy</translation>
    </message>
    <message>
        <source>Overview</source>
        <translation>Prehľad</translation>
    </message>
    <message>
        <source>Legend</source>
        <translation>Legenda</translation>
    </message>
    <message>
        <source>You are using QGIS version %1 built against code revision %2.</source>
        <translation type="unfinished">Používate verziu %1 programu QGIS postavenú na kóde rev. %2.</translation>
    </message>
    <message>
        <source> This copy of QGIS has been built with PostgreSQL support.</source>
        <translation type="unfinished"> Táto kópia QGISu bola postavená s podporou PostgreSQL.</translation>
    </message>
    <message>
        <source> This copy of QGIS has been built without PostgreSQL support.</source>
        <translation type="unfinished"> Táto kópia QGISu bola postavená bez podpory PostgreSQL.</translation>
    </message>
    <message>
        <source>
This binary was compiled against Qt %1,and is currently running against Qt %2</source>
        <translation type="unfinished">
Táto binárna verzia bola skompilovaná s knižnicou Qt %1 a momentálne beží na knižnici Qt %2</translation>
    </message>
    <message>
        <source>This release candidate includes over 120 bug fixes and enchancements over the QGIS 0.9.1 release. In addition we have added the following new features:</source>
        <translation type="obsolete">Tento kandidát na vydanie obsahuje 120 opráv chýb a rozšírení cez vydanie QGIS 0.9.1. Naviac boli pridané nasleduvné funkcie:</translation>
    </message>
    <message>
        <source>Updated icons for improved visual consistancy.</source>
        <translation type="obsolete">Aktualizované ikony pre jednotnejší vzhľad.</translation>
    </message>
    <message>
        <source>Support for migration of old projects to work in newer QGIS versions.</source>
        <translation type="obsolete">Podpora pre prechod starších súborov projektu, tak aby fungovali aj v novších verziách QGISu.</translation>
    </message>
    <message>
        <source>Custom CRS...</source>
        <translation>Vlastný CRS...</translation>
    </message>
    <message>
        <source>Manage custom coordinate reference systems</source>
        <translation type="unfinished">Správa vlastných referenčných súradnicových systémov</translation>
    </message>
    <message>
        <source>Progress bar that displays the status of rendering layers and other time-intensive operations</source>
        <translation type="unfinished">Ukazovateľ priebehu zobrazuje stav vykresľovania vrstiev a iné oparácie citlivé na čas</translation>
    </message>
    <message>
        <source>Shows the map coordinates at the current cursor position. The display is continuously updated as the mouse is moved.</source>
        <translation type="unfinished">Zobrazuje mapové súradnice aktuálnej polohy kurzora. Zobrazenie je pri pohybe myšiu priebežne aktualizované.</translation>
    </message>
    <message>
        <source>Stop map rendering</source>
        <translation>Zastaviť vykresľovanie mapy</translation>
    </message>
    <message>
        <source>This icon shows whether on the fly coordinate reference system transformation is enabled or not. Click the icon to bring up the project properties dialog to alter this behaviour.</source>
        <translation type="unfinished">Táto ikona ukazuje, či je povolený (zapnutý) priamy prevod medzi referenčnými súradnicovými systémami. Kliknutím na túto ikonu sa otvorí dialóg Vlastnosti projektu, kde je možné zmeniť toto nastavenie.</translation>
    </message>
    <message>
        <source>CRS status - Click to open coordinate reference system dialog</source>
        <translation type="unfinished">Stav CRS - kliknutím sa otvorí dialógové okno s nastaveniami referenčého súradnicového systému</translation>
    </message>
    <message>
        <source>Map canvas. This is where raster and vector layers are displayed when added to the map</source>
        <translation type="unfinished">Mapové plátno - je to miesto kde sa zobrazujú rastrové a vektorové vrstvy pridané do mapy</translation>
    </message>
    <message>
        <source>This release candidate includes over 60 bug fixes and enchancements over the QGIS 0.10.0 release. In addition we have added the following new features:</source>
        <translation type="obsolete">Tento kandidát na vydanie obsahuje 60 opráv chýb a rozšírení cez vydanie QGIS 0.10.0. Naviac boli pridané nasledovné funkcie: </translation>
    </message>
    <message>
        <source>Revision of all dialogs for user interface consistancy</source>
        <translation type="obsolete">Revízia všetkých dialógových okien užívateľského rozhrania a ich zjednotenie</translation>
    </message>
    <message>
        <source>Improvements to unique value renderer vector dialog</source>
        <translation type="obsolete">Vylepšenia dialógového okna vykresľovača vektorovej symboliky</translation>
    </message>
    <message>
        <source>Symbol previews when defining vector classes</source>
        <translation type="obsolete">Náhľady na symboly pri definícii tried vo vektorových vrstvách</translation>
    </message>
    <message>
        <source>Separation of python support into its own library</source>
        <translation type="obsolete">Oddelenie podpory Pythonu do samostatnej knižnice</translation>
    </message>
    <message>
        <source>List view and filter for GRASS toolbox to find tools more quickly</source>
        <translation type="obsolete">Typ pohľadu &quot;zoznam&quot; a filtrovanie v Nástrojoch GRASSu, aby bolo možné rýchlešie nájsť požadovaný nástroj</translation>
    </message>
    <message>
        <source>List view and filter for Plugin Manager to find plugins more easily</source>
        <translation type="obsolete">Typ pohľadu &quot;zoznam&quot; a filtrovanie v Správcovi zásuvných modulov, aby bolo možné rýchlešie nájsť požadovaný zásuvný modul</translation>
    </message>
    <message>
        <source>Updated Spatial Reference System definitions</source>
        <translation type="obsolete">Aktualizované definície o priestorových referenčných systémoch</translation>
    </message>
    <message>
        <source>QML Style support for rasters and database layers</source>
        <translation type="obsolete">Podpora štýlov QML pre rastrové a databázové vrstvy</translation>
    </message>
    <message>
        <source>There was an error loading a plugin.The following diagnostic information may help the QGIS developers resolve the issue:
%1.</source>
        <translation type="obsolete">Pri nahrávaní zásuvného modulu sa vyskytla chyba. Nasledujúce diagnostické informácie môžu pomôcť vývojárom QGISu vyriešiť problém:
%1.</translation>
    </message>
    <message>
        <source>Maptips require an active layer</source>
        <translation type="unfinished">Mapové tipy vyžadujú nejakú aktívnu vrstvu</translation>
    </message>
    <message>
        <source>&lt;tt&gt;Settings:Options:General&lt;/tt&gt;</source>
        <comment>Menu path to setting options</comment>
        <translation>&lt;tt&gt;Nastavenia:Vlastnosti:Všeobecné&lt;/tt&gt;</translation>
    </message>
    <message>
        <source>Multiple Instances of QgisApp</source>
        <translation>Viaceré inštancie QgisApp</translation>
    </message>
    <message>
        <source>Multiple instances of Quantum GIS application object detected.
Please contact the developers.
</source>
        <translation>Boli zistené viaceré inštancie objektu aplikácie Quantum GIS.
Prosím kontaktujte vývojárov.
</translation>
    </message>
    <message>
        <source>Shift+Ctrl+S</source>
        <comment>Save Project under a new name</comment>
        <translation type="obsolete">Shift+Ctrl+S</translation>
    </message>
    <message>
        <source>&amp;Print Composer</source>
        <translation>&amp;Skladateľ tlačových výstupov</translation>
    </message>
    <message>
        <source>Ctrl+P</source>
        <comment>Print Composer</comment>
        <translation type="obsolete">Ctrl+P</translation>
    </message>
    <message>
        <source>Print Composer</source>
        <translation>Skladateľ tlačových výstupov</translation>
    </message>
    <message>
        <source>&amp;Undo</source>
        <translation>&amp;Vrátiť späť</translation>
    </message>
    <message>
        <source>Ctrl+Z</source>
        <translation type="unfinished">Ctrl+Z</translation>
    </message>
    <message>
        <source>Undo the last operation</source>
        <translation type="unfinished">Vrátiť späť poslednú operáciu</translation>
    </message>
    <message>
        <source>Cu&amp;t</source>
        <translation>Vys&amp;trihnúť</translation>
    </message>
    <message>
        <source>Ctrl+X</source>
        <translation type="unfinished">Ctrl+X</translation>
    </message>
    <message>
        <source>Cut the current selection&apos;s contents to the clipboard</source>
        <translation type="unfinished">Vystrihne aktuálny obsah a vloží ho do schránky</translation>
    </message>
    <message>
        <source>&amp;Copy</source>
        <translation>&amp;Kopírovať</translation>
    </message>
    <message>
        <source>Ctrl+C</source>
        <translation>Ctrl+C</translation>
    </message>
    <message>
        <source>Copy the current selection&apos;s contents to the clipboard</source>
        <translation type="unfinished">Kopírovať aktuálny obsah výberu do schránky</translation>
    </message>
    <message>
        <source>&amp;Paste</source>
        <translation>&amp;Vložiť</translation>
    </message>
    <message>
        <source>Ctrl+V</source>
        <translation>Ctrl+V</translation>
    </message>
    <message>
        <source>Paste the clipboard&apos;s contents into the current selection</source>
        <translation type="unfinished">Vložiť obsah schránku do aktuálneho výberu</translation>
    </message>
    <message>
        <source>M</source>
        <comment>Measure a Line</comment>
        <translation type="obsolete">M</translation>
    </message>
    <message>
        <source>J</source>
        <comment>Measure an Area</comment>
        <translation type="obsolete">J</translation>
    </message>
    <message>
        <source>Zoom to Selection</source>
        <translation type="unfinished">Na veľkosť výberu</translation>
    </message>
    <message>
        <source>Ctrl+J</source>
        <comment>Zoom to Selection</comment>
        <translation type="obsolete">Ctrl+J</translation>
    </message>
    <message>
        <source>Zoom Actual Size</source>
        <translation type="unfinished">Priblížiť na aktuálnu veľkosť</translation>
    </message>
    <message>
        <source>Zoom to Actual Size</source>
        <translation type="unfinished">Priblížiť na aktuálu veľkosť</translation>
    </message>
    <message>
        <source>Add Vector Layer...</source>
        <translation>Pridať vektorovú vrstvu...</translation>
    </message>
    <message>
        <source>Add Raster Layer...</source>
        <translation>Pridať rastrovú vrstvu...</translation>
    </message>
    <message>
        <source>Add PostGIS Layer...</source>
        <translation>Pridať vrstvu PostGIS...</translation>
    </message>
    <message>
        <source>W</source>
        <comment>Add a Web Mapping Server Layer</comment>
        <translation type="obsolete">W</translation>
    </message>
    <message>
        <source>Add a Web Mapping Server Layer</source>
        <translation type="unfinished">Pridať vrstvu Webovej Mapovej Služby</translation>
    </message>
    <message>
        <source>Open Attribute Table</source>
        <translation>Otvoriť tabuľku atribútov</translation>
    </message>
    <message>
        <source>Save as Shapefile...</source>
        <translation>Uložiť ako súbor shape...</translation>
    </message>
    <message>
        <source>Save the current layer as a shapefile</source>
        <translation type="unfinished">Uloží aktuálnu vrstvu ako súbor shape</translation>
    </message>
    <message>
        <source>Save Selection as Shapefile...</source>
        <translation>Uložiť výber ako súbor shape...</translation>
    </message>
    <message>
        <source>Save the selection as a shapefile</source>
        <translation type="unfinished">Uloží výber ako súbor shape</translation>
    </message>
    <message>
        <source>Properties...</source>
        <translation>Vlastnosti...</translation>
    </message>
    <message>
        <source>Set properties of the current layer</source>
        <translation>Nastaviť vlastnosti aktuálnej vrstvy</translation>
    </message>
    <message>
        <source>Add to Overview</source>
        <translation type="unfinished">Pridať do Prehľadu</translation>
    </message>
    <message>
        <source>Add All to Overview</source>
        <translation type="unfinished">Všetko pridať do Prehľadu</translation>
    </message>
    <message>
        <source>Manage Plugins...</source>
        <translation>Správa zásuvných modulov...</translation>
    </message>
    <message>
        <source>Toggle Full Screen Mode</source>
        <translation>Prepnúť do celoobrazovkového režimu</translation>
    </message>
    <message>
        <source>Minimize</source>
        <translation>Minimalizovať</translation>
    </message>
    <message>
        <source>Ctrl+M</source>
        <comment>Minimize Window</comment>
        <translation type="obsolete">Ctrl+M</translation>
    </message>
    <message>
        <source>Minimizes the active window to the dock</source>
        <translation>Minimalizuje aktívne okno do doku</translation>
    </message>
    <message>
        <source>Zoom</source>
        <translation>Približovanie/oddaľovanie</translation>
    </message>
    <message>
        <source>Toggles between a predefined size and the window size set by the user</source>
        <translation type="unfinished">Prepnúť medzi prednastavenou veľkosťou a veľkosťou okna nastavenou používateľom</translation>
    </message>
    <message>
        <source>Bring All to Front</source>
        <translation>Preniesť všetko do popredia</translation>
    </message>
    <message>
        <source>Bring forward all open windows</source>
        <translation>Preniesť všetky otvorené okná dozadu</translation>
    </message>
    <message>
        <source>&amp;Edit</source>
        <translation>&amp;Upraviť</translation>
    </message>
    <message>
        <source>Panels</source>
        <translation>Panely</translation>
    </message>
    <message>
        <source>Toolbars</source>
        <translation>Panely nástrojov</translation>
    </message>
    <message>
        <source>&amp;Window</source>
        <translation>&amp;Okno</translation>
    </message>
    <message>
        <source>Toggle extents and mouse position display</source>
        <translation type="unfinished">Prepínanie medzi zobrazovaním rozsahu a polohy myši</translation>
    </message>
    <message>
        <source>Choose a file name to save the QGIS project file as</source>
        <translation type="unfinished">Vyberte meno súboru do ktorého sa uloží súbor projektu QGIS</translation>
    </message>
    <message>
        <source>Choose a file name to save the map image as</source>
        <translation type="unfinished">Vyberte meno súboru, do ktorého sa má uložiť obrázok mapy</translation>
    </message>
    <message>
        <source>Start editing failed</source>
        <translation type="unfinished">Pokus o úpravy zlyhal</translation>
    </message>
    <message>
        <source>Provider cannot be opened for editing</source>
        <translation type="unfinished">Správca sa nedá otvoriť pre zápis</translation>
    </message>
    <message>
        <source>Stop editing</source>
        <translation>Ukončiť úpravy</translation>
    </message>
    <message>
        <source>Do you want to save the changes to layer %1?</source>
        <translation>Prajete si uložiť zmeny vo vrstve %1?</translation>
    </message>
    <message>
        <source>Could not commit changes to layer %1

Errors:  %2
</source>
        <translation>Nemožno uložiť zmeny vo vrstve %1

Chyby:  %2
</translation>
    </message>
    <message>
        <source>Problems during roll back</source>
        <translation>Problémy v priebehu návratu do východzieho stavu (roll back)</translation>
    </message>
    <message>
        <source>Python Console</source>
        <translation>Konzola Pythonu</translation>
    </message>
    <message>
        <source>Map coordinates for the current view extents</source>
        <translation type="unfinished">Mapové súradnice pre rozsah aktuálneho pohľadu</translation>
    </message>
    <message>
        <source>This release candidate includes over 265 bug fixes and enchancements over the QGIS 0.11.0 release. In addition we have added the following new features:</source>
        <translation type="obsolete">Tento kandidát na vydanie obsahuje 60 opráv chýb a rozšírení cez vydanie QGIS 0.10.0. Naviac boli pridané nasledovné funkcie:  {265 ?} {0.11.0 ?}</translation>
    </message>
    <message>
        <source>HIG Compliance improvements for Windows / Mac OS X / KDE / Gnome</source>
        <translation type="unfinished">Vylepšené užívateľské rozhranie pre Windows/Mac OS X/KDE/Gnome</translation>
    </message>
    <message>
        <source>Saving a vector layer or subset of that layer to disk with a different Coordinate Reference System to the original.</source>
        <translation type="unfinished">Možnosť uloženia vektorovej vrstvy alebo jej časti na disk v inom než pôdovnom súradnicovom systéme.</translation>
    </message>
    <message>
        <source>Advanced topological editing of vector data.</source>
        <translation type="unfinished">Pokročilé topologické editovanie vektorových údajov.</translation>
    </message>
    <message>
        <source>Single click selection of vector features.</source>
        <translation type="unfinished">Výber vektorových objektov jedným klikom.</translation>
    </message>
    <message>
        <source>Many improvements to raster rendering and support for building pyramids external to the raster file.</source>
        <translation type="unfinished">Viacero vylepšení vykresľovania rastrov a podpora tvorby pyramíd mimo rastrovho súboru.</translation>
    </message>
    <message>
        <source>Overhaul of the map composer for much improved printing support.</source>
        <translation type="unfinished">Oprava Skladateľa mapových výstupov pre ešte lepšiu podporu tlače.</translation>
    </message>
    <message>
        <source>A new &apos;coordinate capture&apos; plugin was added that lets you click on the map and then cut &amp; paste the coordinates to and from the clipboard.</source>
        <translation type="unfinished">Pridaný nový zásuvný modul &apos;Odchyt súradníc&apos; umožňujúci kliknúť do mapy a získať súradnicu, ktorá je zapísaná do schránky a je možné ju získať cez Vystrihnúť &amp; Vložiť.</translation>
    </message>
    <message>
        <source>A new plugin for converting between OGR supported formats was added.</source>
        <translation type="unfinished">Pridaný nový zásuvný modul na prevody medzi formátmi podporovanými knižnicou OGR.</translation>
    </message>
    <message>
        <source>A new plugin for converting from DXF files to shapefiles was added.</source>
        <translation type="unfinished">Pridaný nový zásuvný modul na prevody súborov DXF do formátu Shape.</translation>
    </message>
    <message>
        <source>A new plugin was added for interpolating point features into ASCII grid layers.</source>
        <translation type="unfinished">Pridaný nový zásuvný modul na interpoláciu bodových objektov do vrstiev ASCII grid.</translation>
    </message>
    <message>
        <source>Plugin toolbar positions are now correctly saved when the application is closed.</source>
        <translation type="unfinished">Poloha panelov s nástrojmi je odteraz pri zavretí aplikácie ukladaná správne.</translation>
    </message>
    <message>
        <source>In the WMS client, WMS standards support has been improved.</source>
        <translation type="unfinished">V klientovi WMS bola vylepšená podpora štandardu WMS.</translation>
    </message>
    <message>
        <source>Complete API revision - we now have a stable API following well defined naming conventions.</source>
        <translation>Úplná revízia API - stabilné API s dobre definovanými mennými konvenciami.</translation>
    </message>
    <message>
        <source>Ported all GDAL/OGR and GEOS usage to use C APIs only.</source>
        <translation type="unfinished">Portované použitie knižníc GDAL/OGR a GEOS len do C-čkového API.</translation>
    </message>
    <message>
        <source>The python plugin installer was completely overhauled, the new version having many improvements, including checking that the version of QGIS running will support a plugin that is being installed.</source>
        <translation type="unfinished">Inštalátor zásuvných modulov v jazyku Python bolo úplne prepracované - nová verzia má mnoho vylepšení vrátane kontroly, či je inštalovaný modul vhodný pre práve bežiacu verziu QGISu.</translation>
    </message>
    <message>
        <source>Vector editing overhaul - handling of geometry and attribute edit transactions is now handled transparently in one place.</source>
        <translation type="unfinished">Prepracovanie úprav vektorových vrstiev - obsluha geometrie a úravy atribútov sú prehľadne obsluhované z jedného miesta.</translation>
    </message>
    <message>
        <source>Changes</source>
        <translation>Zmeny</translation>
    </message>
    <message>
        <source>QGIS 1.0.2 is a bug fix release for the stable version of QGIS.A summary of the improvements can be found at https://trac.osgeo.org/qgis/query?status=closed&amp;milestone=Version+1.0.2</source>
        <translation>QGIS 1.0.2 je opravné vydanie stabilnej verzie QGISu. Súhrn vylepšení možno nájsť na stránke https://trac.osgeo.org/qgis/query?status=closed&amp;milestone=Version+1.0.2</translation>
    </message>
    <message>
        <source>QGIS 1.0.1 is a bug fix release for the stable version of QGIS.A summary of the improvements can be found at https://trac.osgeo.org/qgis/query?status=closed&amp;milestone=Version+1.0.1</source>
        <translation>QGIS 1.0.1 je opravné vydanie stabilnej verzie QGISu. Súhrn vylepšení možno nájsť na stránke https://trac.osgeo.org/qgis/query?status=closed&amp;milestone=Version+1.0.1</translation>
    </message>
    <message>
        <source>The QGIS 1.0 release includes over 265 bug fixes and enchancements over the QGIS 0.11.0 release. In addition we have added the following new features:</source>
        <translation type="obsolete">Vydanie QGISu verzie 1.0 obsahuje viac ako 265 opráv a rozšírení oproti verzii QGIS 0.11.0. Naviac boli pridané nasledujúce nové vlastnosti:</translation>
    </message>
    <message>
        <source>The QGIS 1.0 release includes over 265 bug fixes and enhancements over the QGIS 0.11.0 release. In addition we have added the following new features:</source>
        <translation type="unfinished"></translation>
    </message>
    <message>
        <source>Ctrl+N</source>
        <translation type="unfinished">Ctrl+N</translation>
    </message>
    <message>
        <source>Ctrl+O</source>
        <translation type="unfinished">Ctrl+O</translation>
    </message>
    <message>
        <source>Ctrl+S</source>
        <translation type="unfinished">Ctrl+S</translation>
    </message>
    <message>
        <source>Shift+Ctrl+S</source>
        <translation type="unfinished">Shift+Ctrl+S</translation>
    </message>
    <message>
        <source>Ctrl+I</source>
        <translation type="unfinished">Ctrl+I</translation>
    </message>
    <message>
        <source>Ctrl+P</source>
        <translation type="unfinished">Ctrl+P</translation>
    </message>
    <message>
        <source>Ctrl+Q</source>
        <translation type="unfinished">Ctrl+Q</translation>
    </message>
    <message>
        <source>.</source>
        <translation type="unfinished">.</translation>
    </message>
    <message>
        <source>/</source>
        <translation type="unfinished">/</translation>
    </message>
    <message>
        <source>Ctrl+/</source>
        <translation type="unfinished">Ctrl+/</translation>
    </message>
    <message>
        <source>Ctrl++</source>
        <translation type="unfinished">Ctrl++</translation>
    </message>
    <message>
        <source>Ctrl+-</source>
        <translation type="unfinished">Ctrl+-</translation>
    </message>
    <message>
        <source>I</source>
        <translation type="unfinished">I</translation>
    </message>
    <message>
        <source>M</source>
        <translation type="unfinished">M</translation>
    </message>
    <message>
        <source>J</source>
        <translation type="unfinished">J</translation>
    </message>
    <message>
        <source>F</source>
        <translation type="unfinished">F</translation>
    </message>
    <message>
        <source>Ctrl+J</source>
        <translation type="unfinished">Ctrl+J</translation>
    </message>
    <message>
        <source>Ctrl+B</source>
        <translation type="unfinished">Ctrl+B</translation>
    </message>
    <message>
        <source>B</source>
        <translation type="unfinished">B</translation>
    </message>
    <message>
        <source>Ctrl+R</source>
        <translation type="unfinished">Ctrl+R</translation>
    </message>
    <message>
        <source>N</source>
        <translation type="unfinished"></translation>
    </message>
    <message>
        <source>V</source>
        <translation type="unfinished">V</translation>
    </message>
    <message>
        <source>R</source>
        <translation type="unfinished">R</translation>
    </message>
    <message>
        <source>D</source>
        <translation type="unfinished">D</translation>
    </message>
    <message>
        <source>W</source>
        <translation type="unfinished"></translation>
    </message>
    <message>
        <source>Ctrl+D</source>
        <translation type="unfinished">Ctrl+D</translation>
    </message>
    <message>
        <source>O</source>
        <translation type="unfinished">O</translation>
    </message>
    <message>
        <source>+</source>
        <translation type="unfinished">+</translation>
    </message>
    <message>
        <source>-</source>
        <translation type="unfinished">-</translation>
    </message>
    <message>
        <source>S</source>
        <translation type="unfinished"></translation>
    </message>
    <message>
        <source>H</source>
        <translation type="unfinished">H</translation>
    </message>
    <message>
        <source>Ctrl-F</source>
        <translation type="unfinished">Ctrl-F</translation>
    </message>
    <message>
        <source>P</source>
        <translation type="unfinished">P</translation>
    </message>
    <message>
        <source>Ctrl+M</source>
        <translation type="unfinished">Ctrl+M</translation>
    </message>
    <message>
        <source>Ctrl+?</source>
        <translation type="unfinished">Ctrl+?</translation>
    </message>
    <message>
        <source>F1</source>
        <translation type="unfinished">F1</translation>
    </message>
    <message>
        <source>Ctrl+H</source>
        <translation type="unfinished">Ctrl+H</translation>
    </message>
    <message>
        <source></source>
        <translation></translation>
    </message>
</context>
<context>
    <name>QgisAppBase</name>
    <message>
        <source>MainWindow</source>
        <translation type="obsolete">HlavneOkno</translation>
    </message>
    <message>
        <source>Legend</source>
        <translation type="obsolete">Legenda</translation>
    </message>
    <message>
        <source>Map View</source>
        <translation type="obsolete">Mapový pohľad</translation>
    </message>
    <message>
        <source>QGIS</source>
        <translation>QGIS</translation>
    </message>
</context>
<context>
    <name>QgsAbout</name>
    <message>
        <source>About Quantum GIS</source>
        <translation>O programe Quantum GIS</translation>
    </message>
    <message>
        <source>Ok</source>
        <translation>OK</translation>
    </message>
    <message>
        <source>About</source>
        <translation>O programe</translation>
    </message>
    <message>
        <source>Version</source>
        <translation>Verzia</translation>
    </message>
    <message>
        <source>What&apos;s New</source>
        <translation>Čo je nového</translation>
    </message>
    <message>
        <source>http://www.gnu.org/licenses</source>
        <translation>http://www.gnu.org/licenses</translation>
    </message>
    <message>
        <source>QGIS Home Page</source>
        <translation>Domovská stránka QGIS</translation>
    </message>
    <message>
        <source>Providers</source>
        <translation>Prístup k údajom</translation>
    </message>
    <message>
        <source>Developers</source>
        <translation>Vývojári</translation>
    </message>
    <message>
        <source>Name</source>
        <translation type="obsolete">Meno</translation>
    </message>
    <message>
        <source>Website</source>
        <translation type="obsolete">Domovská stránka</translation>
    </message>
    <message>
        <source>Sponsors</source>
        <translation>Sponzori</translation>
    </message>
    <message>
        <source>Quantum GIS is licensed under the GNU General Public License</source>
        <translation>Quantum GIS je šírený pod licenciou GNU General Public License</translation>
    </message>
    <message>
        <source>&lt;p&gt;The following have sponsored QGIS by contributing money to fund development and other project costs&lt;/p&gt;</source>
        <translation type="obsolete">&lt;p&gt;Nasledujúci ľudia podporili QGIS financovaním jeho vývoja a ďalších nákladov projektu&lt;/p&gt;</translation>
    </message>
    <message>
        <source>Available QGIS Data Provider Plugins</source>
        <translation type="obsolete">Dostupné zásuvné moduly QGISu na spracovanie údajov</translation>
    </message>
    <message>
        <source>Available Qt Database Plugins</source>
        <translation type="obsolete">Dostupné zásuvné moduly Qt na prácu s databázou</translation>
    </message>
    <message>
        <source>Available Qt Image Plugins</source>
        <translation type="obsolete">Dostupné zásuvné moduly Qt na prácu s obrazom</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Lucida Grande&apos;; font-size:13pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:16px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:x-large; font-weight:600;&quot;&gt;&lt;span style=&quot; font-size:x-large;&quot;&gt;Quantum GIS (QGIS)&lt;/span&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Lucida Grande&apos;; font-size:13pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:16px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:x-large; font-weight:600;&quot;&gt;&lt;span style=&quot; font-size:x-large;&quot;&gt;Quantum GIS (QGIS)&lt;/span&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>Join our user mailing list</source>
        <translation>Pripojte sa k našej e-mailovej konferencii</translation>
    </message>
</context>
<context>
    <name>QgsAddAttrDialogBase</name>
    <message>
        <source>Add Attribute</source>
        <translation>Pridať atribút</translation>
    </message>
    <message>
        <source>Name:</source>
        <translation>Meno:</translation>
    </message>
    <message>
        <source>Type:</source>
        <translation>Typ:</translation>
    </message>
</context>
<context>
    <name>QgsApplication</name>
    <message>
        <source>Exception</source>
        <translation>Výnimka</translation>
    </message>
</context>
<context>
    <name>QgsAttributeActionDialog</name>
    <message>
        <source>Select an action</source>
        <comment>File dialog window title</comment>
        <translation>Výber akcie</translation>
    </message>
</context>
<context>
    <name>QgsAttributeActionDialogBase</name>
    <message>
        <source>Name</source>
        <translation>názov</translation>
    </message>
    <message>
        <source>Action</source>
        <translation>Akcia</translation>
    </message>
    <message>
        <source>Move up</source>
        <translation>Posunúť nahor</translation>
    </message>
    <message>
        <source>Move the selected action up</source>
        <translation>Posunie vybratú akciu nahor</translation>
    </message>
    <message>
        <source>Move down</source>
        <translation>Posunúť nadol</translation>
    </message>
    <message>
        <source>Move the selected action down</source>
        <translation>Posunie vybratú akciu nadol</translation>
    </message>
    <message>
        <source>Remove</source>
        <translation>Odobrať</translation>
    </message>
    <message>
        <source>Remove the selected action</source>
        <translation>Odoberie vybratú akciu</translation>
    </message>
    <message>
        <source>Enter the name of an action here. The name should be unique (qgis will make it unique if necessary).</source>
        <translation>Sem vložte meno akcie. Názov by mal byť jedinečný (ak to bude potrebné, qgis ho upraví, aby bol jedinečný).</translation>
    </message>
    <message>
        <source>Enter the action name here</source>
        <translation>Sem vložte názov akcie</translation>
    </message>
    <message>
        <source>Enter the action command here</source>
        <translation>Sem vložte príkaz akcie</translation>
    </message>
    <message>
        <source>Insert action</source>
        <translation>Vložiť akciu</translation>
    </message>
    <message>
        <source>Inserts the action into the list above</source>
        <translation>Do zoznamu hore vloží akciu</translation>
    </message>
    <message>
        <source>Update action</source>
        <translation>Aktualizovať akciu</translation>
    </message>
    <message>
        <source>Update the selected action</source>
        <translation>Aktualizuje vybranú akciu</translation>
    </message>
    <message>
        <source>Insert field</source>
        <translation>Vložiť pole</translation>
    </message>
    <message>
        <source>Inserts the selected field into the action, prepended with a %</source>
        <translation>Vloží vybrané pole do akcie, predchádzané s %</translation>
    </message>
    <message>
        <source>The valid attribute names for this layer</source>
        <translation>Platné mená atribútov pre túto vrstvu</translation>
    </message>
    <message>
        <source>This list contains all actions that have been defined for the current layer. Add actions by entering the details in the controls below and then pressing the Insert action button. Actions can be edited here by double clicking on the item.</source>
        <translation>Tento zoznam obsahuje všetky akcie, ktoré boli definované pre aktuálnu vrstvu. Akcie možno pridávať vyplnením príslušných políčok a následným stlačením tlačidla Vložiť akciu. Akcie je možné upraviť dvojklikom na príslušnú položku.</translation>
    </message>
    <message>
        <source>Capture</source>
        <translation>Sledovanie</translation>
    </message>
    <message>
        <source>Capture output</source>
        <translation>Sledovať výstup</translation>
    </message>
    <message>
        <source>Captures any output from the action</source>
        <translation>Sleduje akýkoľvek výstup z akcie</translation>
    </message>
    <message>
        <source>Captures the standard output or error generated by the action and displays it in a dialog box</source>
        <translation>Sleduje štandardný výstup alebo chybu generovanú akciou a zobrazí ju v dialógovom okne</translation>
    </message>
    <message>
        <source>Enter the action here. This can be any program, script or command that is available on your system. When the action is invoked any set of characters that start with a % and then have the name of a field will be replaced by the value of that field. The special characters %% will be replaced by the value of the field that was selected. Double quote marks group text into single arguments to the program, script or command. Double quotes will be ignored if preceeded by a backslash</source>
        <translation>Sem vložte akciu. Môže ňou byť ľubovoľný program, skript, alebo príkaz, ktorý je dostupný vo vašom systéme. Keď je vyvolaná akcia, akákoľvek skupina znakov začínajúca so znakom % nasledovaným menom poľa bude nahradená hodnotou poľa. Špeciálne znaky %% budú nahradené hodnotou poľa, ktoré bolo vybraté. Dvojité úvodzovky zoskupia text do jednotlivých argumentov programu, skriptu alebo príkazu. Dvojité úvodzovky budú ignorované pokiaľ je pred nimi spätná lomka</translation>
    </message>
    <message>
        <source>Attribute Actions</source>
        <translation>Akcie atribútov</translation>
    </message>
    <message>
        <source>Action properties</source>
        <translation>Vlastnosti akcie</translation>
    </message>
    <message>
        <source>Browse for action</source>
        <translation>Vybrať akciu</translation>
    </message>
    <message>
        <source>Click to browse for an action</source>
        <translation>Kliknutím môžete vyhľadať akciu</translation>
    </message>
    <message>
        <source>Clicking the buttone will let you select an application to use as the action</source>
        <translation type="obsolete">Po kliknutí na tlačidlo je možné vybrať aplikáciu ktorá sa bude používaž pre túto akciu</translation>
    </message>
    <message>
        <source>...</source>
        <translation>...</translation>
    </message>
    <message>
        <source>Clicking the button will let you select an application to use as the action</source>
        <translation>Po kliknutí na tlačidlo je možné vybrať aplikáciu, ktorá sa použije ako akcia</translation>
    </message>
</context>
<context>
    <name>QgsAttributeDialog</name>
    <message>
        <source> (int)</source>
        <translation> (int)</translation>
    </message>
    <message>
        <source> (dbl)</source>
        <translation> (dbl)</translation>
    </message>
    <message>
        <source> (txt)</source>
        <translation> (txt)</translation>
    </message>
    <message>
        <source>...</source>
        <translation>...</translation>
    </message>
    <message>
        <source>Select a file</source>
        <translation>Výber súboru</translation>
    </message>
</context>
<context>
    <name>QgsAttributeDialogBase</name>
    <message>
        <source>Enter Attribute Values</source>
        <translation>Vložiť hodnoty atribútov</translation>
    </message>
    <message>
        <source>1</source>
        <translation type="obsolete">1</translation>
    </message>
    <message>
        <source>Attribute</source>
        <translation type="obsolete">Atribút</translation>
    </message>
    <message>
        <source>Value</source>
        <translation type="obsolete">Hodnota</translation>
    </message>
</context>
<context>
    <name>QgsAttributeTable</name>
    <message>
        <source>Run action</source>
        <translation>Spustiť akciu</translation>
    </message>
    <message>
        <source>Updating selection...</source>
        <translation>Aktualizuje sa výber...</translation>
    </message>
    <message>
        <source>Abort</source>
        <translation>Prerušiť</translation>
    </message>
</context>
<context>
    <name>QgsAttributeTableBase</name>
    <message>
        <source>Attribute Table</source>
        <translation>Tabuľka atribútov</translation>
    </message>
    <message>
        <source>Start editing</source>
        <translation type="obsolete">Začať upravovať</translation>
    </message>
    <message>
        <source>&amp;Close</source>
        <translation type="obsolete">&amp;Zatvoriť</translation>
    </message>
    <message>
        <source>Alt+C</source>
        <translation type="obsolete">Alt+Z</translation>
    </message>
    <message>
        <source>Ctrl+X</source>
        <translation type="obsolete">Ctrl+X</translation>
    </message>
    <message>
        <source>Ctrl+N</source>
        <translation type="obsolete">Ctrl+N</translation>
    </message>
    <message>
        <source>Ctrl+S</source>
        <translation>Ctrl+S</translation>
    </message>
    <message>
        <source>Invert selection</source>
        <translation>Invertovať výber</translation>
    </message>
    <message>
        <source>Ctrl+T</source>
        <translation>Ctrl+T</translation>
    </message>
    <message>
        <source>Move selected to top</source>
        <translation>Presunúť výber navrch</translation>
    </message>
    <message>
        <source>Remove selection</source>
        <translation>Odobrať výber</translation>
    </message>
    <message>
        <source>Copy selected rows to clipboard (Ctrl+C)</source>
        <translation>Skopírovať vybrané riadky do schránky (Ctrl+C)</translation>
    </message>
    <message>
        <source>Copies the selected rows to the clipboard</source>
        <translation>Skopíruje vybrané riadky do schránky</translation>
    </message>
    <message>
        <source>Ctrl+C</source>
        <translation>Ctrl+C</translation>
    </message>
    <message>
        <source>Stop editin&amp;g</source>
        <translation type="obsolete">U&amp;končiť úpravy</translation>
    </message>
    <message>
        <source>Alt+G</source>
        <translation type="obsolete">Alt+K</translation>
    </message>
    <message>
        <source>Search for:</source>
        <translation type="obsolete">Nájsť: </translation>
    </message>
    <message>
        <source>in</source>
        <translation>v</translation>
    </message>
    <message>
        <source>Search</source>
        <translation>Nájsť</translation>
    </message>
    <message>
        <source>Adva&amp;nced...</source>
        <translation>&amp;Pokročilé...</translation>
    </message>
    <message>
        <source>Alt+N</source>
        <translation>Alt+N</translation>
    </message>
    <message>
        <source>&amp;Help</source>
        <translation type="obsolete">&amp;Pomocník</translation>
    </message>
    <message>
        <source>New column</source>
        <translation type="obsolete">Nový stĺpec</translation>
    </message>
    <message>
        <source>Delete column</source>
        <translation type="obsolete">Vymazať stĺpec</translation>
    </message>
    <message>
        <source>Zoom map to the selected rows (Ctrl-F)</source>
        <translation type="obsolete">Zmeniť pohľad mapy na veľkosť vybraných riadkov (Ctrl+F)</translation>
    </message>
    <message>
        <source>Zoom map to the selected rows</source>
        <translation type="unfinished">Nastaví mapový pohľad na veľkosť objektov zodpovedajúcich vybraným riadkom</translation>
    </message>
    <message>
        <source>Ctrl+F</source>
        <translation type="obsolete">Ctrl+F</translation>
    </message>
    <message>
        <source>Toggle editing mode</source>
        <translation>Prepnutie na režim úprav</translation>
    </message>
    <message>
        <source>Click to toggle table editing</source>
        <translation>Po kliknutí sa tabuľka prepne do režimu úprav</translation>
    </message>
    <message>
        <source>Search for</source>
        <translation>Hľadať</translation>
    </message>
    <message>
        <source>Zoom map to the selected rows (Ctrl-J)</source>
        <translation>Zmeniť pohľad mapy na veľkosť vybraných riadkov (Ctrl+J)</translation>
    </message>
    <message>
        <source>Ctrl+J</source>
        <translation>Ctrl+J</translation>
    </message>
</context>
<context>
    <name>QgsAttributeTableDisplay</name>
    <message>
        <source>select</source>
        <translation>označiť</translation>
    </message>
    <message>
        <source>select and bring to top</source>
        <translation>označiť a preniesť navrch</translation>
    </message>
    <message>
        <source>show only matching</source>
        <translation>ukázať len zodpovedajúce</translation>
    </message>
    <message>
        <source>Search string parsing error</source>
        <translation>Chyba pri kontrole reťazca</translation>
    </message>
    <message>
        <source>Search results</source>
        <translation>Výsledky hľadania</translation>
    </message>
    <message>
        <source>You&apos;ve supplied an empty search string.</source>
        <translation>Zadali ste prázdny reťazec.</translation>
    </message>
    <message>
        <source>Error during search</source>
        <translation>Chyba počas hľadania</translation>
    </message>
    <message>
        <source>Found %d matching features.</source>
        <translation type="obsolete">
        </translation>
    </message>
    <message>
        <source>No matching features found.</source>
        <translation>Nenašli sa žiadne zodpovedajúce objekty.</translation>
    </message>
    <message>
        <source>Name conflict</source>
        <translation type="obsolete">Konflikt názov</translation>
    </message>
    <message>
        <source>Stop editing</source>
        <translation type="obsolete">Ukončiť úpravy</translation>
    </message>
    <message>
        <source>Do you want to save the changes?</source>
        <translation type="obsolete">Prajete si uložiť zmeny?</translation>
    </message>
    <message>
        <source>Error</source>
        <translation type="obsolete">Chyba</translation>
    </message>
    <message>
        <source>The attribute could not be inserted. The name already exists in the table.</source>
        <translation type="obsolete">Atribút nemožno vložiť. Atribút s takýmto názvom už v tabuľke existuje.</translation>
    </message>
    <message>
        <source>Could not commit changes - changes are still pending</source>
        <translation type="obsolete">Nemožno odoslať zmeny - zmeny sa stále vykonávajú</translation>
    </message>
    <message>
        <source>Editing not permitted</source>
        <translation type="obsolete">Úpravy nie sú povolené</translation>
    </message>
    <message>
        <source>The data provider is read only, editing is not allowed.</source>
        <translation type="obsolete">Správca údajov umožňuje len čítanie, úpravy nie sú povolené.</translation>
    </message>
    <message>
        <source>Start editing</source>
        <translation type="obsolete">Začať upravovať</translation>
    </message>
    <message>
        <source>Attribute table - </source>
        <translation>Tabuľka atribútov -</translation>
    </message>
    <message>
        <source>QGIS</source>
        <translation>QGIS</translation>
    </message>
    <message>
        <source>File</source>
        <translation>Súbor</translation>
    </message>
    <message>
        <source>Close</source>
        <translation>Zatvoriť</translation>
    </message>
    <message>
        <source>Ctrl+W</source>
        <translation>Ctrl+W</translation>
    </message>
    <message>
        <source>Edit</source>
        <translation>Upraviť</translation>
    </message>
    <message>
        <source>&amp;Undo</source>
        <translation>&amp;Vrátiť späť</translation>
    </message>
    <message>
        <source>Ctrl+Z</source>
        <translation type="unfinished">Ctrl+Z</translation>
    </message>
    <message>
        <source>Cu&amp;t</source>
        <translation>Vys&amp;trihnúť</translation>
    </message>
    <message>
        <source>Ctrl+X</source>
        <translation>Ctrl+X</translation>
    </message>
    <message>
        <source>&amp;Copy</source>
        <translation>&amp;Kopírovať</translation>
    </message>
    <message>
        <source>Ctrl+C</source>
        <translation>Ctrl+C</translation>
    </message>
    <message>
        <source>&amp;Paste</source>
        <translation>&amp;Vložiť</translation>
    </message>
    <message>
        <source>Ctrl+V</source>
        <translation>Ctrl+V</translation>
    </message>
    <message>
        <source>Delete</source>
        <translation>Vymazať</translation>
    </message>
    <message>
        <source>Layer</source>
        <translation>Vrstva</translation>
    </message>
    <message>
        <source>Zoom to Selection</source>
        <translation type="unfinished">Na veľkosť výberu</translation>
    </message>
    <message>
        <source>Ctrl+J</source>
        <translation type="unfinished">Ctrl+J</translation>
    </message>
    <message>
        <source>Toggle Editing</source>
        <translation>Prepnúť na úpravy</translation>
    </message>
    <message>
        <source>Table</source>
        <translation>Tabuľka</translation>
    </message>
    <message>
        <source>Move to Top</source>
        <translation>Presunúť navrch</translation>
    </message>
    <message>
        <source>Invert</source>
        <translation>Invertovať</translation>
    </message>
    <message>
        <source>bad_alloc exception</source>
        <translation>výnimka bad_alloc</translation>
    </message>
    <message>
        <source>Filling the attribute table has been stopped because there was no more virtual memory left</source>
        <translation>Naplnenie tabuľky atribútov bolo zastavené, pretože nebol dostatok virtuálnej pamäte</translation>
    </message>
    <message>
        <source>Found %1 matching features.</source>
        <translation type="obsolete">Našlo sa %1 zodpovedajúcich objektov.
        </translation>
    </message>
</context>
<context>
    <name>QgsBookmarks</name>
    <message>
        <source>Really Delete?</source>
        <translation>Skutočne zmazať?</translation>
    </message>
    <message>
        <source>Are you sure you want to delete the </source>
        <translation>Ste si istý, že chcete vymazať záložku </translation>
    </message>
    <message>
        <source> bookmark?</source>
        <translation>?</translation>
    </message>
    <message>
        <source>Error deleting bookmark</source>
        <translation>Chyba pri odstraňovaní záložky</translation>
    </message>
    <message>
        <source>Failed to delete the </source>
        <translation> Zlyhal pokus o odstránenie záložky</translation>
    </message>
    <message>
        <source> bookmark from the database. The database said:
</source>
        <translation> z databázy. Správa z databázy:
</translation>
    </message>
    <message>
        <source>&amp;Delete</source>
        <translation>&amp;Vymazať</translation>
    </message>
    <message>
        <source>&amp;Zoom to</source>
        <translation>&amp;Priblížiť na</translation>
    </message>
</context>
<context>
    <name>QgsBookmarksBase</name>
    <message>
        <source>Geospatial Bookmarks</source>
        <translation>Geopriestorové záložky</translation>
    </message>
    <message>
        <source>Name</source>
        <translation>Meno</translation>
    </message>
    <message>
        <source>Project</source>
        <translation>Projekt</translation>
    </message>
    <message>
        <source>Extent</source>
        <translation>Rozsah</translation>
    </message>
    <message>
        <source>Id</source>
        <translation>Id</translation>
    </message>
</context>
<context>
    <name>QgsComposer</name>
    <message>
        <source>Choose a filename to save the map image as</source>
        <translation type="obsolete">Vyberte meno súboru, do ktorého sa má uložiť obrázok mapy</translation>
    </message>
    <message>
        <source>Choose a filename to save the map as</source>
        <translation type="obsolete">Vyberte meno súboru do ktorého sa má uložiť mapa</translation>
    </message>
    <message>
        <source> for read/write</source>
        <translation type="obsolete">na čítanie/zápis</translation>
    </message>
    <message>
        <source>Error in Print</source>
        <translation type="obsolete">Chyba pri tlači</translation>
    </message>
    <message>
        <source>Cannot seek</source>
        <translation type="obsolete">Nemožno požadovať</translation>
    </message>
    <message>
        <source>Cannot overwrite BoundingBox</source>
        <translation type="obsolete">Nemožno prepísať BoundingBox (obdĺžnik ohraničenia)</translation>
    </message>
    <message>
        <source>Cannot find BoundingBox</source>
        <translation type="obsolete">Nemožno nájsť BoundingBox (obdĺžnik ohraničenia)</translation>
    </message>
    <message>
        <source>Cannot overwrite translate</source>
        <translation type="obsolete">Nemožno prepísať prevod</translation>
    </message>
    <message>
        <source>Cannot find translate</source>
        <translation type="obsolete">Nemožno nájsť prevod</translation>
    </message>
    <message>
        <source>File IO Error</source>
        <translation type="obsolete">Vstupno/výstupná chyba súboru</translation>
    </message>
    <message>
        <source>Paper does not match</source>
        <translation type="obsolete">Nezodpovedá veľkosť papiera</translation>
    </message>
    <message>
        <source>The selected paper size does not match the composition size</source>
        <translation type="obsolete">Vybraná veľkosť papiera nezodpovedá veľkosti kompozície</translation>
    </message>
    <message>
        <source>Big image</source>
        <translation>Priveľký obrázok</translation>
    </message>
    <message>
        <source>To create image </source>
        <translation> Na vytvorenie obrázka</translation>
    </message>
    <message>
        <source> requires circa </source>
        <translation>  je potrebných približne</translation>
    </message>
    <message>
        <source> MB of memory</source>
        <translation> MB pamäte</translation>
    </message>
    <message>
        <source>SVG warning</source>
        <translation>Upozornenie (SVG)</translation>
    </message>
    <message>
        <source>Don&apos;t show this message again</source>
        <translation>Nabudúce už túto správu nezobrazovať</translation>
    </message>
    <message>
        <source>QGIS - print composer</source>
        <translation>QGIS - Skladateľ tlačových výstupov</translation>
    </message>
    <message>
        <source>Map 1</source>
        <translation>Mapa 1</translation>
    </message>
    <message>
        <source>Couldn&apos;t open </source>
        <translation type="obsolete">Nemožno otvoriť </translation>
    </message>
    <message>
        <source>format</source>
        <translation>formát</translation>
    </message>
    <message>
        <source>SVG Format</source>
        <translation>Formát SVG</translation>
    </message>
    <message>
        <source>&lt;p&gt;The SVG export function in Qgis has several problems due to bugs and deficiencies in the </source>
        <translation>&lt;p&gt;Funkcia exportu do SVG je v QGISe dosť probémová kvôli chybám a nedostatkom v </translation>
    </message>
    <message>
        <source>Move Content</source>
        <translation type="obsolete">Presunúť obsah</translation>
    </message>
    <message>
        <source>Move item content</source>
        <translation type="obsolete">Presunúť obsah položiek</translation>
    </message>
    <message>
        <source>&amp;Group</source>
        <translation type="obsolete">&amp;Zoskupiť</translation>
    </message>
    <message>
        <source>Group items</source>
        <translation type="obsolete">Zoskupiť položky</translation>
    </message>
    <message>
        <source>&amp;Ungroup</source>
        <translation type="obsolete">Zr&amp;ušiť zoskupenie</translation>
    </message>
    <message>
        <source>Ungroup items</source>
        <translation type="obsolete">Zrušiť zoskupenie položiek</translation>
    </message>
    <message>
        <source>Raise</source>
        <translation type="obsolete">Vyššie</translation>
    </message>
    <message>
        <source>Raise selected items</source>
        <translation type="obsolete">Posunúť vybrané položky vyššie</translation>
    </message>
    <message>
        <source>Lower</source>
        <translation type="obsolete">Nižšie</translation>
    </message>
    <message>
        <source>Lower selected items</source>
        <translation type="obsolete">Posunúť vybrané položky nižšie</translation>
    </message>
    <message>
        <source>Bring to Front</source>
        <translation type="obsolete">Preniesť dopredu</translation>
    </message>
    <message>
        <source>Move selected items to top</source>
        <translation type="obsolete">Preniesť vybrané položky navrch</translation>
    </message>
    <message>
        <source>Send to Back</source>
        <translation type="obsolete">Presniesť dozadu</translation>
    </message>
    <message>
        <source>Move selected items to bottom</source>
        <translation type="obsolete">Preniesť vybrané položky naspodok</translation>
    </message>
    <message>
        <source>QGIS</source>
        <translation>QGIS</translation>
    </message>
    <message>
        <source>File</source>
        <translation>Súbor</translation>
    </message>
    <message>
        <source>Close</source>
        <translation>Zatvoriť</translation>
    </message>
    <message>
        <source>Ctrl+W</source>
        <translation>Ctrl+W</translation>
    </message>
    <message>
        <source>Edit</source>
        <translation>Upraviť</translation>
    </message>
    <message>
        <source>&amp;Undo</source>
        <translation>&amp;Vrátiť späť</translation>
    </message>
    <message>
        <source>Ctrl+Z</source>
        <translation>Ctrl+Z</translation>
    </message>
    <message>
        <source>Cu&amp;t</source>
        <translation>Vys&amp;trihnúť</translation>
    </message>
    <message>
        <source>Ctrl+X</source>
        <translation>Ctrl+X</translation>
    </message>
    <message>
        <source>&amp;Copy</source>
        <translation>&amp;Kopírovať</translation>
    </message>
    <message>
        <source>Ctrl+C</source>
        <translation>Ctrl+C</translation>
    </message>
    <message>
        <source>&amp;Paste</source>
        <translation>&amp;Vložiť</translation>
    </message>
    <message>
        <source>Ctrl+V</source>
        <translation>Ctrl+V</translation>
    </message>
    <message>
        <source>Delete</source>
        <translation>Vymazať</translation>
    </message>
    <message>
        <source>View</source>
        <translation>Pohľad</translation>
    </message>
    <message>
        <source>Layout</source>
        <translation>Rozvrhnutie</translation>
    </message>
    <message>
        <source>Choose a file name to save the map image as</source>
        <translation type="unfinished">Vyberte meno súboru, do ktorého sa má uložiť obrázok mapy</translation>
    </message>
    <message>
        <source>Choose a file name to save the map as</source>
        <translation type="unfinished">Vyberte meno súboru do ktorého sa má uložiť mapa</translation>
    </message>
    <message>
        <source>Project contains WMS layers</source>
        <translation>Projekt obsahuje vrstvy WMS </translation>
    </message>
    <message>
        <source>Some WMS servers (e.g. UMN mapserver) have a limit for the WIDTH and HEIGHT parameter. Printing layers from such servers may exceed this limit. If this is the case, the WMS layer will not be printed</source>
        <translation>Niektoré WMS servery (napr. UMN mapserver) majú obmedzený parameter WIDTH a HEIGHT. Vrstvy na tlačenie z takýchto serverov môžu presiahnuť toto obmedzenie. V tomto prípade nebude WMS vrstva vytlačená</translation>
    </message>
    <message>
        <source>Qt4 svg code. Of note, text does not appear in the SVG file and there are problems with the map bounding box clipping other items such as the legend or scale bar.&lt;/p&gt;</source>
        <translation>kóde knižnice Qt4 pre prácu s svg. Napríklad v súbore SVG neobjaví text a vyskytujú sa aj problémy s obdĺžnikom záujmového územia obsahujúcim iné položky ako je legenda alebo grafická mierka.&lt;/p&gt;</translation>
    </message>
    <message>
        <source>Qt4 svg code. In particular, there are problems with layers not being clipped to the map bounding box.&lt;/p&gt;</source>
        <translation>kóde knižnice Qt4 pre prácu s svg. Konkréne ide o problémy s vrstvami ktoré nie sú orezané podľa obdĺžnika záujmového územia.&lt;/p&gt;</translation>
    </message>
    <message>
        <source>If you require a vector-based output file from Qgis it is suggested that you try printing to PostScript if the SVG output is not satisfactory.&lt;/p&gt;</source>
        <translation>Pokiaľ potrebujete výstupný súbor založený na vekotorvej grafike a výstup vo fromáte SVG nie je uspokojivý, odporúča sa tlač do formátu PostScript.&lt;/p&gt;</translation>
    </message>
</context>
<context>
    <name>QgsComposerBase</name>
    <message>
        <source>General</source>
        <translation>Všeobecné</translation>
    </message>
    <message>
        <source>Composition</source>
        <translation>Kompozícia</translation>
    </message>
    <message>
        <source>Item</source>
        <translation>Položka</translation>
    </message>
    <message>
        <source>&amp;Open Template ...</source>
        <translation type="obsolete">&amp;Otvoriť šablónu ...</translation>
    </message>
    <message>
        <source>Save Template &amp;As...</source>
        <translation type="obsolete">Uložiť šablónu &amp;ako...</translation>
    </message>
    <message>
        <source>&amp;Print...</source>
        <translation>&amp;Tlačiť...</translation>
    </message>
    <message>
        <source>Add new map</source>
        <translation>Pridá novú mapu</translation>
    </message>
    <message>
        <source>Add new label</source>
        <translation>Pridá nový popis</translation>
    </message>
    <message>
        <source>Add new vect legend</source>
        <translation>Pridá novú vektorovú legendu</translation>
    </message>
    <message>
        <source>Select/Move item</source>
        <translation>Vybrať/premiestniť položku</translation>
    </message>
    <message>
        <source>Export as image</source>
        <translation type="obsolete">Exportovať ako obrázok</translation>
    </message>
    <message>
        <source>Export as SVG</source>
        <translation type="obsolete">Export do formátu SVG</translation>
    </message>
    <message>
        <source>Add new scalebar</source>
        <translation>Pridať novú grafickú mierku</translation>
    </message>
    <message>
        <source>Refresh view</source>
        <translation>Obnoviť pohľad</translation>
    </message>
    <message>
        <source>MainWindow</source>
        <translation>Skladateľ máp</translation>
    </message>
    <message>
        <source>Zoom All</source>
        <translation type="obsolete">Celá strana</translation>
    </message>
    <message>
        <source>Zoom In</source>
        <translation>Priblížiť</translation>
    </message>
    <message>
        <source>Zoom Out</source>
        <translation>Oddialiť</translation>
    </message>
    <message>
        <source>Add Image</source>
        <translation>Pridať obrázok</translation>
    </message>
    <message>
        <source>Close</source>
        <translation>Zatvoriť</translation>
    </message>
    <message>
        <source>Help</source>
        <translation>Pomocník</translation>
    </message>
    <message>
        <source>&amp;Open Template...</source>
        <translation type="obsolete">&amp;Otvoriť šablónu...</translation>
    </message>
    <message>
        <source>Zoom Full</source>
        <translation>Celá mapa</translation>
    </message>
    <message>
        <source>Add Map</source>
        <translation>Pridať mapu</translation>
    </message>
    <message>
        <source>Add Label</source>
        <translation>Pridať popis</translation>
    </message>
    <message>
        <source>Add Vector Legend</source>
        <translation>Pridať vektorovú legendu</translation>
    </message>
    <message>
        <source>Move Item</source>
        <translation>Posunúť položku</translation>
    </message>
    <message>
        <source>Export as Image...</source>
        <translation>Exportovať ako obrázok...</translation>
    </message>
    <message>
        <source>Export as SVG...</source>
        <translation>Export do formátu SVG...</translation>
    </message>
    <message>
        <source>Add Scalebar</source>
        <translation>Pridať grafickú mierku</translation>
    </message>
    <message>
        <source>Refresh</source>
        <translation>Obnoviť</translation>
    </message>
    <message>
        <source>Ungroup</source>
        <translation>Zrušiť zoskupenie</translation>
    </message>
    <message>
        <source>Raise</source>
        <translation>Vyššie</translation>
    </message>
    <message>
        <source>Lower</source>
        <translation>Nižšie</translation>
    </message>
    <message>
        <source>Move Content</source>
        <translation>Presunúť obsah</translation>
    </message>
    <message>
        <source>Move item content</source>
        <translation>Presunúť obsah položiek</translation>
    </message>
    <message>
        <source>Group</source>
        <translation>Zoskupiť</translation>
    </message>
    <message>
        <source>Group items</source>
        <translation>Zoskupiť položky</translation>
    </message>
    <message>
        <source>Ungroup items</source>
        <translation>Zrušiť zoskupenie položiek</translation>
    </message>
    <message>
        <source>Raise selected items</source>
        <translation>Posunúť vybrané položky vyššie</translation>
    </message>
    <message>
        <source>Lower selected items</source>
        <translation>Posunúť vybrané položky nižšie</translation>
    </message>
    <message>
        <source>Bring to Front</source>
        <translation>Preniesť dopredu</translation>
    </message>
    <message>
        <source>Move selected items to top</source>
        <translation>Preniesť vybrané položky navrch</translation>
    </message>
    <message>
        <source>Send to Back</source>
        <translation>Preniesť dozadu</translation>
    </message>
    <message>
        <source>Move selected items to bottom</source>
        <translation>Preniesť vybrané položky naspodok</translation>
    </message>
</context>
<context>
    <name>QgsComposerItemWidgetBase</name>
    <message>
        <source>Form</source>
        <translation>Vlastnosti mapovej položky</translation>
    </message>
    <message>
        <source>Composer item properties</source>
        <translation>Vlastnosti položky</translation>
    </message>
    <message>
        <source>Color:</source>
        <translation>Farba:</translation>
    </message>
    <message>
        <source>Frame...</source>
        <translation>Farba rámu...</translation>
    </message>
    <message>
        <source>Background...</source>
        <translation>Farba pozadia...</translation>
    </message>
    <message>
        <source>Opacity:</source>
        <translation>Viditeľnosť:</translation>
    </message>
    <message>
        <source>Outline width: </source>
        <translation>Hrúbka obrysu: </translation>
    </message>
    <message>
        <source>Frame</source>
        <translation>Rám</translation>
    </message>
</context>
<context>
    <name>QgsComposerLabelBase</name>
    <message>
        <source>Label Options</source>
        <translation type="obsolete">Nastavenia popisu</translation>
    </message>
    <message>
        <source>Font</source>
        <translation type="obsolete">Písmo</translation>
    </message>
    <message>
        <source>Box</source>
        <translation type="obsolete">Rám</translation>
    </message>
</context>
<context>
    <name>QgsComposerLabelWidgetBase</name>
    <message>
        <source>Label Options</source>
        <translation>Vlastnosti popisu</translation>
    </message>
    <message>
        <source>Font</source>
        <translation>Písmo</translation>
    </message>
    <message>
        <source>Margin (mm):</source>
        <translation>Okraj (mm):</translation>
    </message>
</context>
<context>
    <name>QgsComposerLegendItemDialogBase</name>
    <message>
        <source>Legend item properties</source>
        <translation>Vlastnosti položky legendy</translation>
    </message>
    <message>
        <source>Item text:</source>
        <translation>Text položky:</translation>
    </message>
</context>
<context>
    <name>QgsComposerLegendWidgetBase</name>
    <message>
        <source>Barscale Options</source>
        <translation>Vlastnosti grafickej mierky</translation>
    </message>
    <message>
        <source>General</source>
        <translation>Všeobecné</translation>
    </message>
    <message>
        <source>Title:</source>
        <translation>Titulok:</translation>
    </message>
    <message>
        <source>Font:</source>
        <translation>Písmo:</translation>
    </message>
    <message>
        <source>Title...</source>
        <translation>Titulok...</translation>
    </message>
    <message>
        <source>Layer...</source>
        <translation>Vrstva...</translation>
    </message>
    <message>
        <source>Item...</source>
        <translation>Položka...</translation>
    </message>
    <message>
        <source>Symbol width: </source>
        <translation>Šírka symbolu: </translation>
    </message>
    <message>
        <source>Symbol height:</source>
        <translation>Výška symbolu:</translation>
    </message>
    <message>
        <source>Layer space: </source>
        <translation>Priestor vrstvy: </translation>
    </message>
    <message>
        <source>Symbol space:</source>
        <translation>Priestor pre symbol:</translation>
    </message>
    <message>
        <source>Icon label space:</source>
        <translation>Priestor pre popis:</translation>
    </message>
    <message>
        <source>Box space:</source>
        <translation type="unfinished">Priestor pre rám:</translation>
    </message>
    <message>
        <source>Legend items</source>
        <translation type="unfinished">Položky legendy</translation>
    </message>
    <message>
        <source>down</source>
        <translation>nižšie</translation>
    </message>
    <message>
        <source>up</source>
        <translation>vyššie</translation>
    </message>
    <message>
        <source>remove</source>
        <translation>odobrať</translation>
    </message>
    <message>
        <source>edit...</source>
        <translation>upraviť...</translation>
    </message>
    <message>
        <source>update</source>
        <translation>aktualizovať</translation>
    </message>
    <message>
        <source>update all</source>
        <translation>aktualizovať všetky</translation>
    </message>
</context>
<context>
    <name>QgsComposerMap</name>
    <message>
        <source>Extent (calculate scale)</source>
        <translation type="obsolete">Rozsah (vypočíta mierku)</translation>
    </message>
    <message>
        <source>Scale (calculate extent)</source>
        <translation type="obsolete">Mierka (vypočíta rozsah)</translation>
    </message>
    <message>
        <source>Map %1</source>
        <translation type="obsolete">Mapa %1</translation>
    </message>
    <message>
        <source>Cache</source>
        <translation type="obsolete">Z vyrovnávacej pamäte</translation>
    </message>
    <message>
        <source>Render</source>
        <translation type="obsolete">Prekresľovať</translation>
    </message>
    <message>
        <source>Rectangle</source>
        <translation type="obsolete">Iba obdĺžnik</translation>
    </message>
    <message>
        <source>Map</source>
        <translation>Mapa</translation>
    </message>
    <message>
        <source>Map will be printed here</source>
        <translation>Tu bude vytlačená mapa</translation>
    </message>
</context>
<context>
    <name>QgsComposerMapBase</name>
    <message>
        <source>Map options</source>
        <translation type="obsolete">Nastavenia mapy</translation>
    </message>
    <message>
        <source>&lt;b&gt;Map&lt;/b&gt;</source>
        <translation type="obsolete">&lt;b&gt;Mapa&lt;/b&gt;</translation>
    </message>
    <message>
        <source>Width</source>
        <translation type="obsolete">Šírka</translation>
    </message>
    <message>
        <source>Height</source>
        <translation type="obsolete">Výška</translation>
    </message>
    <message>
        <source>Set Extent</source>
        <translation type="obsolete">Nastaviť rozsah</translation>
    </message>
    <message>
        <source>Set map extent to current extent in QGIS map canvas</source>
        <translation type="obsolete">Nastaví rozsah mapy podľa aktuálneho rozsahu na mapovom plátne QGISu</translation>
    </message>
    <message>
        <source>Line width scale</source>
        <translation type="obsolete">Mierka hrúbky čiary</translation>
    </message>
    <message>
        <source>Width of one unit in millimeters</source>
        <translation type="obsolete">Jednotkou hrúbky je milimeter</translation>
    </message>
    <message>
        <source>Symbol scale</source>
        <translation type="obsolete">Mierka symbolov</translation>
    </message>
    <message>
        <source>Font size scale</source>
        <translation type="obsolete">Mierka veľkosti písma</translation>
    </message>
    <message>
        <source>Frame</source>
        <translation type="obsolete">Rám</translation>
    </message>
    <message>
        <source>Preview</source>
        <translation type="obsolete">Náhľad</translation>
    </message>
    <message>
        <source>Set</source>
        <translation type="obsolete">Nastaviť</translation>
    </message>
    <message>
        <source>1:</source>
        <translation type="obsolete">1:</translation>
    </message>
    <message>
        <source>Scale:</source>
        <translation type="obsolete">Mierka:</translation>
    </message>
</context>
<context>
    <name>QgsComposerMapWidget</name>
    <message>
        <source>Cache</source>
        <translation>Z vyrovnávacej pamäte</translation>
    </message>
    <message>
        <source>Rectangle</source>
        <translation>Iba obdĺžnik</translation>
    </message>
    <message>
        <source>Render</source>
        <translation>Prekresľovať</translation>
    </message>
</context>
<context>
    <name>QgsComposerMapWidgetBase</name>
    <message>
        <source>Map options</source>
        <translation>Vlastnosi mapy</translation>
    </message>
    <message>
        <source>&lt;b&gt;Map&lt;/b&gt;</source>
        <translation>&lt;b&gt;Mapa&lt;/b&gt;</translation>
    </message>
    <message>
        <source>Width</source>
        <translation>Šírka</translation>
    </message>
    <message>
        <source>Height</source>
        <translation>Výška</translation>
    </message>
    <message>
        <source>Scale:</source>
        <translation>Mierka:</translation>
    </message>
    <message>
        <source>1:</source>
        <translation>1:</translation>
    </message>
    <message>
        <source>Map extent</source>
        <translation>Rozsah mapy</translation>
    </message>
    <message>
        <source>X min:</source>
        <translation>X min:</translation>
    </message>
    <message>
        <source>Y min:</source>
        <translation>Y min:</translation>
    </message>
    <message>
        <source>X max:</source>
        <translation>X max:</translation>
    </message>
    <message>
        <source>Y max:</source>
        <translation>Y max:</translation>
    </message>
    <message>
        <source>set to map canvas extent</source>
        <translation type="unfinished">nastaviť na rozsah mapového plátna</translation>
    </message>
    <message>
        <source>Preview</source>
        <translation>Náhľad</translation>
    </message>
    <message>
        <source>Update preview</source>
        <translation>Aktualizovať náhľad</translation>
    </message>
</context>
<context>
    <name>QgsComposerPicture</name>
    <message>
        <source>Warning</source>
        <translation type="obsolete">Upozornenie</translation>
    </message>
    <message>
        <source>Cannot load picture.</source>
        <translation type="obsolete">Nemožno nahrať obrázok.</translation>
    </message>
    <message>
        <source>Choose a file</source>
        <translation type="obsolete">Vybrať súbor</translation>
    </message>
    <message>
        <source>Pictures (</source>
        <translation type="obsolete">Obrázky (</translation>
    </message>
</context>
<context>
    <name>QgsComposerPictureBase</name>
    <message>
        <source>Picture Options</source>
        <translation type="obsolete">Vlastnosti obrázka</translation>
    </message>
    <message>
        <source>Frame</source>
        <translation type="obsolete">Rám</translation>
    </message>
    <message>
        <source>Angle</source>
        <translation type="obsolete">Uhol</translation>
    </message>
    <message>
        <source>Width</source>
        <translation type="obsolete">Šírka</translation>
    </message>
    <message>
        <source>Height</source>
        <translation type="obsolete">Výška</translation>
    </message>
    <message>
        <source>Browse</source>
        <translation type="obsolete">Prechádzať</translation>
    </message>
</context>
<context>
    <name>QgsComposerPictureWidget</name>
    <message>
        <source>Select svg or image file</source>
        <translation type="unfinished">Vyberte súbor svg alebo obrázok</translation>
    </message>
</context>
<context>
    <name>QgsComposerPictureWidgetBase</name>
    <message>
        <source>Picture Options</source>
        <translation>Vlastnosti obrázka</translation>
    </message>
    <message>
        <source>Browse...</source>
        <translation>Prechádzať...</translation>
    </message>
    <message>
        <source>Width:</source>
        <translation>Šírka:</translation>
    </message>
    <message>
        <source>Height:</source>
        <translation>Výška:</translation>
    </message>
    <message>
        <source>Rotation:</source>
        <translation>Otočenie:</translation>
    </message>
</context>
<context>
    <name>QgsComposerScaleBar</name>
    <message>
        <source>Single Box</source>
        <translation type="obsolete">Jednoduchý rám</translation>
    </message>
    <message>
        <source>Double Box</source>
        <translation type="obsolete">Dvojitý rám</translation>
    </message>
    <message>
        <source>Line Ticks Middle</source>
        <translation type="obsolete">Čiarky v strede</translation>
    </message>
    <message>
        <source>Line Ticks Down</source>
        <translation type="obsolete">Čiarky dole</translation>
    </message>
    <message>
        <source>Line Ticks Up</source>
        <translation type="obsolete">Čiarky hore</translation>
    </message>
    <message>
        <source>Numeric</source>
        <translation type="obsolete">Číselná</translation>
    </message>
</context>
<context>
    <name>QgsComposerScaleBarWidget</name>
    <message>
        <source>Single Box</source>
        <translation>Jednoduchý obdĺžnik</translation>
    </message>
    <message>
        <source>Double Box</source>
        <translation>Dvojitý obdĺžnik</translation>
    </message>
    <message>
        <source>Line Ticks Middle</source>
        <translation>Čiarky v strede</translation>
    </message>
    <message>
        <source>Line Ticks Down</source>
        <translation>Čiarky dole</translation>
    </message>
    <message>
        <source>Line Ticks Up</source>
        <translation>Čiarky hore</translation>
    </message>
    <message>
        <source>Numeric</source>
        <translation>Číselná</translation>
    </message>
    <message>
        <source>Map </source>
        <translation type="obsolete">Mapa </translation>
    </message>
    <message>
        <source>Map %1</source>
        <translation>Mapa %1</translation>
    </message>
</context>
<context>
    <name>QgsComposerScaleBarWidgetBase</name>
    <message>
        <source>Barscale Options</source>
        <translation>Vlastnosti grafickej mierky</translation>
    </message>
    <message>
        <source>Segment size (map units):</source>
        <translation>Veľkosť úseku (v mapových jednotkách):</translation>
    </message>
    <message>
        <source>Map units per bar unit:</source>
        <translation type="unfinished">Poč. map. jedn. na jedn. graf. mierky:</translation>
    </message>
    <message>
        <source>Number of segments:</source>
        <translation>Počet úsekov:</translation>
    </message>
    <message>
        <source>Segments left:</source>
        <translation>Záporné úseky:</translation>
    </message>
    <message>
        <source>Style:</source>
        <translation>Štýl:</translation>
    </message>
    <message>
        <source>Map:</source>
        <translation>Mapa:</translation>
    </message>
    <message>
        <source>Height (mm):</source>
        <translation>Výška (mm):</translation>
    </message>
    <message>
        <source>Line width:</source>
        <translation>Hrúbka čiary:</translation>
    </message>
    <message>
        <source>Label space:</source>
        <translation>Priestor okolo popisu:</translation>
    </message>
    <message>
        <source>Box space:</source>
        <translation>Priestor rámu:</translation>
    </message>
    <message>
        <source>Unit label:</source>
        <translation type="unfinished">Jednotka popisu:</translation>
    </message>
    <message>
        <source>Font...</source>
        <translation>Písmo...</translation>
    </message>
    <message>
        <source>Color...</source>
        <translation>Farba...</translation>
    </message>
</context>
<context>
    <name>QgsComposerScalebarBase</name>
    <message>
        <source>Barscale Options</source>
        <translation type="obsolete">Možnosti grafickej mierky</translation>
    </message>
    <message>
        <source>Segment size</source>
        <translation type="obsolete">Veľkosť dielika</translation>
    </message>
    <message>
        <source>Number of segments</source>
        <translation type="obsolete">Počet dielikov</translation>
    </message>
    <message>
        <source>Map units per scalebar unit</source>
        <translation type="obsolete">Počet map. jednot. na jednot. graf. mierky</translation>
    </message>
    <message>
        <source>Unit label</source>
        <translation type="obsolete">Značka jednotky</translation>
    </message>
    <message>
        <source>Map</source>
        <translation type="obsolete">Mapa</translation>
    </message>
    <message>
        <source>Font</source>
        <translation type="obsolete">Písmo</translation>
    </message>
    <message>
        <source>Line width</source>
        <translation type="obsolete">Hrúbka čiary</translation>
    </message>
</context>
<context>
    <name>QgsComposerVectorLegend</name>
    <message>
        <source>Combine selected layers</source>
        <translation type="obsolete">Zlúčiť vybrané vrstvy</translation>
    </message>
    <message>
        <source>Cache</source>
        <translation type="obsolete">Z vyrovnávacej pamäte</translation>
    </message>
    <message>
        <source>Render</source>
        <translation type="obsolete">Prekresľovať</translation>
    </message>
    <message>
        <source>Rectangle</source>
        <translation type="obsolete">Iba obdĺžnik</translation>
    </message>
    <message>
        <source>Legend</source>
        <translation type="obsolete">Legenda</translation>
    </message>
</context>
<context>
    <name>QgsComposerVectorLegendBase</name>
    <message>
        <source>Vector Legend Options</source>
        <translation>Možnosti vektorovej legendy</translation>
    </message>
    <message>
        <source>Title</source>
        <translation>Titulok</translation>
    </message>
    <message>
        <source>Map</source>
        <translation>Mapa</translation>
    </message>
    <message>
        <source>Font</source>
        <translation>Písmo</translation>
    </message>
    <message>
        <source>Preview</source>
        <translation>Náhľad</translation>
    </message>
    <message>
        <source>Box</source>
        <translation>Rám</translation>
    </message>
    <message>
        <source>Layers</source>
        <translation>Vrstvy</translation>
    </message>
    <message>
        <source>Group</source>
        <translation>Skupina</translation>
    </message>
    <message>
        <source>ID</source>
        <translation>ID</translation>
    </message>
</context>
<context>
    <name>QgsComposition</name>
    <message>
        <source>Custom</source>
        <translation type="obsolete">Vlastné</translation>
    </message>
    <message>
        <source>A5 (148x210 mm)</source>
        <translation type="obsolete">A5 (148x210 mm)</translation>
    </message>
    <message>
        <source>A4 (210x297 mm)</source>
        <translation type="obsolete">A4 (210x297 mm)</translation>
    </message>
    <message>
        <source>A3 (297x420 mm)</source>
        <translation type="obsolete">A3 (297x420 mm)</translation>
    </message>
    <message>
        <source>A2 (420x594 mm)</source>
        <translation type="obsolete">A2 (420x594 mm)</translation>
    </message>
    <message>
        <source>A1 (594x841 mm)</source>
        <translation type="obsolete">A1 (594x841 mm)</translation>
    </message>
    <message>
        <source>A0 (841x1189 mm)</source>
        <translation type="obsolete">A0 (841x1189 mm)</translation>
    </message>
    <message>
        <source>B5 (176 x 250 mm)</source>
        <translation type="obsolete">B5 (176 x 250 mm)</translation>
    </message>
    <message>
        <source>B4 (250 x 353 mm)</source>
        <translation type="obsolete">B4 (250 x 353 mm)</translation>
    </message>
    <message>
        <source>B3 (353 x 500 mm)</source>
        <translation type="obsolete">B3 (353 x 500 mm)</translation>
    </message>
    <message>
        <source>B2 (500 x 707 mm)</source>
        <translation type="obsolete">B2 (500 x 707 mm)</translation>
    </message>
    <message>
        <source>B1 (707 x 1000 mm)</source>
        <translation type="obsolete">B1 (707 x 1000 mm)</translation>
    </message>
    <message>
        <source>B0 (1000 x 1414 mm)</source>
        <translation type="obsolete">B0 (1000 x 1414 mm)</translation>
    </message>
    <message>
        <source>Letter (8.5x11 inches)</source>
        <translation type="obsolete">Letter (8.5x11 palcov)</translation>
    </message>
    <message>
        <source>Legal (8.5x14 inches)</source>
        <translation type="obsolete">Legal (8.5x14 palcov)</translation>
    </message>
    <message>
        <source>Portrait</source>
        <translation type="obsolete">Na výšku</translation>
    </message>
    <message>
        <source>Landscape</source>
        <translation type="obsolete">Na šírku</translation>
    </message>
    <message>
        <source>Out of memory</source>
        <translation type="obsolete">Nedostaok voľnej pamäte</translation>
    </message>
    <message>
        <source>Qgis is unable to resize the paper size due to insufficient memory.
 It is best that you avoid using the map composer until you restart qgis.
</source>
        <translation type="obsolete"> Qgis nemôže kvôli nedostatku pamäte zmeniť veľkosť papiera.
 Najlepšie bude vyhnúť sa ďalšej práci so Skladateľom máp, kým nereštartujete Qgis.
</translation>
    </message>
    <message>
        <source>Label</source>
        <translation type="obsolete">Popis</translation>
    </message>
    <message>
        <source>Warning</source>
        <translation type="obsolete">Upozornenie</translation>
    </message>
    <message>
        <source>Cannot load picture.</source>
        <translation type="obsolete">Nemožno nahrať obrázok.</translation>
    </message>
</context>
<context>
    <name>QgsCompositionBase</name>
    <message>
        <source>Composition</source>
        <translation>Kompozícia</translation>
    </message>
    <message>
        <source>Paper</source>
        <translation>Papier</translation>
    </message>
    <message>
        <source>Size</source>
        <translation>Veľkosť</translation>
    </message>
    <message>
        <source>Units</source>
        <translation>Jednotky</translation>
    </message>
    <message>
        <source>Width</source>
        <translation>Šírka</translation>
    </message>
    <message>
        <source>Height</source>
        <translation>Výška</translation>
    </message>
    <message>
        <source>Orientation</source>
        <translation>Orientácia</translation>
    </message>
    <message>
        <source>Resolution (dpi)</source>
        <translation type="obsolete">Rozlíšenie (dpi)</translation>
    </message>
</context>
<context>
    <name>QgsCompositionWidget</name>
    <message>
        <source>Landscape</source>
        <translation>Na šírku</translation>
    </message>
    <message>
        <source>Portrait</source>
        <translation>Na výšku</translation>
    </message>
    <message>
        <source>Custom</source>
        <translation>Vlastné</translation>
    </message>
    <message>
        <source>A5 (148x210 mm)</source>
        <translation>A5 (148x210 mm)</translation>
    </message>
    <message>
        <source>A4 (210x297 mm)</source>
        <translation>A4 (210x297 mm)</translation>
    </message>
    <message>
        <source>A3 (297x420 mm)</source>
        <translation>A3 (297x420 mm)</translation>
    </message>
    <message>
        <source>A2 (420x594 mm)</source>
        <translation>A2 (420x594 mm)</translation>
    </message>
    <message>
        <source>A1 (594x841 mm)</source>
        <translation>A1 (594x841 mm)</translation>
    </message>
    <message>
        <source>A0 (841x1189 mm)</source>
        <translation>A0 (841x1189 mm)</translation>
    </message>
    <message>
        <source>B5 (176 x 250 mm)</source>
        <translation>B5 (176 x 250 mm)</translation>
    </message>
    <message>
        <source>B4 (250 x 353 mm)</source>
        <translation>B4 (250 x 353 mm)</translation>
    </message>
    <message>
        <source>B3 (353 x 500 mm)</source>
        <translation>B3 (353 x 500 mm)</translation>
    </message>
    <message>
        <source>B2 (500 x 707 mm)</source>
        <translation>B2 (500 x 707 mm)</translation>
    </message>
    <message>
        <source>B1 (707 x 1000 mm)</source>
        <translation>B1 (707 x 1000 mm)</translation>
    </message>
    <message>
        <source>B0 (1000 x 1414 mm)</source>
        <translation>B0 (1000 x 1414 mm)</translation>
    </message>
    <message>
        <source>Letter (8.5x11 inches)</source>
        <translation>Letter (8.5x11 palcov)</translation>
    </message>
    <message>
        <source>Legal (8.5x14 inches)</source>
        <translation>Legal (8.5x14 palcov)</translation>
    </message>
</context>
<context>
    <name>QgsCompositionWidgetBase</name>
    <message>
        <source>Composition</source>
        <translation>Zostava</translation>
    </message>
    <message>
        <source>Paper</source>
        <translation>Papier</translation>
    </message>
    <message>
        <source>Orientation</source>
        <translation>Orientácia</translation>
    </message>
    <message>
        <source>Height</source>
        <translation>Výška</translation>
    </message>
    <message>
        <source>Width</source>
        <translation>Šírka</translation>
    </message>
    <message>
        <source>Units</source>
        <translation>Jednotky</translation>
    </message>
    <message>
        <source>Size</source>
        <translation>Veľkosť</translation>
    </message>
    <message>
        <source>Print quality (dpi)</source>
        <translation>Kvalita tlače (dpi)</translation>
    </message>
</context>
<context>
    <name>QgsConnectionDialogBase</name>
    <message>
        <source>Connection Information</source>
        <translation type="obsolete">Informácie o spojení</translation>
    </message>
    <message>
        <source>Host</source>
        <translation type="obsolete">Hostiteľ</translation>
    </message>
    <message>
        <source>Database</source>
        <translation type="obsolete">Databáza</translation>
    </message>
    <message>
        <source>Username</source>
        <translation type="obsolete">Meno používateľa</translation>
    </message>
    <message>
        <source>Name</source>
        <translation type="obsolete">Meno</translation>
    </message>
    <message>
        <source>Name of the new connection</source>
        <translation type="obsolete">Meno nového spojenia</translation>
    </message>
    <message>
        <source>Password</source>
        <translation type="obsolete">Heslo</translation>
    </message>
    <message>
        <source>Test Connect</source>
        <translation type="obsolete">Vyskúšať spojenie</translation>
    </message>
    <message>
        <source>Save Password</source>
        <translation type="obsolete">Uložiť heslo</translation>
    </message>
    <message>
        <source>Create a New PostGIS connection</source>
        <translation type="obsolete">Vytvoriť nové spojenie k PostGIS </translation>
    </message>
    <message>
        <source>Port</source>
        <translation type="obsolete">Port</translation>
    </message>
    <message>
        <source>5432</source>
        <translation type="obsolete">5432</translation>
    </message>
</context>
<context>
    <name>QgsContinuousColorDialogBase</name>
    <message>
        <source>Continuous color</source>
        <translation>Spojitá farba</translation>
    </message>
    <message>
        <source>Maximum Value:</source>
        <translation>Maximálna hodnota:</translation>
    </message>
    <message>
        <source>Outline Width:</source>
        <translation>Hrúbka obrysu:</translation>
    </message>
    <message>
        <source>Minimum Value:</source>
        <translation>Minimálna hodnota:</translation>
    </message>
    <message>
        <source>Classification Field:</source>
        <translation>Triediť podľa poľa:</translation>
    </message>
    <message>
        <source>Draw polygon outline</source>
        <translation>Vykresľovať obrys polygónu</translation>
    </message>
</context>
<context>
    <name>QgsCoordinateTransform</name>
    <message>
        <source>Failed</source>
        <translation>Zlyhala</translation>
    </message>
    <message>
        <source>transform of</source>
        <translation>transformácia</translation>
    </message>
    <message>
        <source>with error: </source>
        <translation>s nasledujúcou chybou:</translation>
    </message>
    <message>
        <source>The source spatial reference system (SRS) is not valid. </source>
        <translation type="obsolete">Neplatný zdrojový priestorový referenčný systém (SRS). </translation>
    </message>
    <message>
        <source>The coordinates can not be reprojected. The SRS is: </source>
        <translation type="obsolete">Súradnice nemožno previesť. SRS je:</translation>
    </message>
    <message>
        <source>The destination spatial reference system (SRS) is not valid. </source>
        <translation type="obsolete">Cieľový priestorový referenčný systém (SRS) nie je platný. </translation>
    </message>
    <message>
        <source>The source spatial reference system (CRS) is not valid. </source>
        <translation type="unfinished">Neplatný zdrojový priestorový referenčný systém (CRS). </translation>
    </message>
    <message>
        <source>The coordinates can not be reprojected. The CRS is: </source>
        <translation type="unfinished">Súradnice nemožno previesť. CRS je: </translation>
    </message>
    <message>
        <source>The destination spatial reference system (CRS) is not valid. </source>
        <translation type="unfinished">Cieľový priestorový referenčný systém (SRS) nie je platný. </translation>
    </message>
</context>
<context>
    <name>QgsCopyrightLabelPlugin</name>
    <message>
        <source>Bottom Left</source>
        <translation>Vľavo dole</translation>
    </message>
    <message>
        <source>Top Left</source>
        <translation>Vľavo hore</translation>
    </message>
    <message>
        <source>Top Right</source>
        <translation>Vpravo hore</translation>
    </message>
    <message>
        <source>&amp;Decorations</source>
        <translation>&amp;Doplnky</translation>
    </message>
    <message>
        <source>Creates a copyright label that is displayed on the map canvas.</source>
        <translation>Vytvorí označenie autorských práv zobrazujúce sa na mapovom plátne.</translation>
    </message>
    <message>
        <source>Bottom Right</source>
        <translation>Vpravo dole</translation>
    </message>
    <message>
        <source>&amp;Copyright Label</source>
        <translation>&amp;Označenie autorských práv</translation>
    </message>
</context>
<context>
    <name>QgsCopyrightLabelPluginGuiBase</name>
    <message>
        <source>Copyright Label Plugin</source>
        <translation>Zásuvný modul na označenie copyrightu</translation>
    </message>
    <message>
        <source>Placement</source>
        <translation>Umiestnenie</translation>
    </message>
    <message>
        <source>Bottom Left</source>
        <translation>Vľavo dole</translation>
    </message>
    <message>
        <source>Top Left</source>
        <translation>Vľavo hore</translation>
    </message>
    <message>
        <source>Bottom Right</source>
        <translation>Vpravo dole</translation>
    </message>
    <message>
        <source>Top Right</source>
        <translation>Vpravo hore</translation>
    </message>
    <message>
        <source>Orientation</source>
        <translation>Orientácia</translation>
    </message>
    <message>
        <source>Horizontal</source>
        <translation>Vodorovne</translation>
    </message>
    <message>
        <source>Vertical</source>
        <translation>Zvisle</translation>
    </message>
    <message>
        <source>Enable Copyright Label</source>
        <translation>Zapnúť označenie copyrightu</translation>
    </message>
    <message>
        <source>Color</source>
        <translation>Farba</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;span style=&quot; font-size:12pt;&quot;&gt;Description&lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot;-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Enter your copyright label below. This plugin supports basic html markup tags for formatting the label. For example:&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;span style=&quot; font-weight:600;&quot;&gt;&amp;lt;B&amp;gt; Bold text &amp;lt;/B&amp;gt; &lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-weight:600;&quot;&gt;&lt;span style=&quot; font-weight:400; font-style:italic;&quot;&gt;&amp;lt;I&amp;gt; Italics &amp;lt;/I&amp;gt;&lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-style:italic;&quot;&gt;&lt;span style=&quot; font-style:normal;&quot;&gt;(note: &amp;amp;copy; gives a copyright symbol)&lt;/span&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;span style=&quot; font-size:12pt;&quot;&gt;Popis&lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot;-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Do poľa dolu vložte text vášho copyrightu (označenia vlastníckych práv). Tento zásuvný modul podporuje pri formátovaní základné značky jazyka HTML. Napríklad:&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;span style=&quot; font-weight:600;&quot;&gt;&amp;lt;B&amp;gt; Tučný text &amp;lt;/B&amp;gt; &lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-weight:600;&quot;&gt;&lt;span style=&quot; font-weight:400; font-style:italic;&quot;&gt;&amp;lt;I&amp;gt; Kurzíva &amp;lt;/I&amp;gt;&lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-style:italic;&quot;&gt;&lt;span style=&quot; font-style:normal;&quot;&gt;(poznámka: &amp;amp;copy; sa používa pre znak copyrightu)&lt;/span&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>&#xa9; QGIS 2008</source>
        <translation type="obsolete">© QGIS 2008</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;&lt;span style=&quot; font-size:12pt;&quot;&gt;Description&lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot;-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;Enter your copyright label below. This plugin supports basic html markup tags for formatting the label. For example:&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;&lt;span style=&quot; font-weight:600;&quot;&gt;&amp;lt;B&amp;gt; Bold text &amp;lt;/B&amp;gt; &lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt; font-weight:600;&quot;&gt;&lt;span style=&quot; font-weight:400; font-style:italic;&quot;&gt;&amp;lt;I&amp;gt; Italics &amp;lt;/I&amp;gt;&lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt; font-style:italic;&quot;&gt;&lt;span style=&quot; font-style:normal;&quot;&gt;(note: &amp;amp;copy; gives a copyright symbol)&lt;/span&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;&lt;span style=&quot; font-size:12pt;&quot;&gt;Popis&lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot;-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;Do poľa pod týmto textom zadajte text vášho copyrightu (označenia vlastníckych práv). Tento zásuvný modul dokáže rozpoznať pri formátovaní základné značky jazyka HTML. Napríklad:&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;&lt;span style=&quot; font-weight:600;&quot;&gt;&amp;lt;B&amp;gt; Tučný text &amp;lt;/B&amp;gt; &lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt; font-weight:600;&quot;&gt;&lt;span style=&quot; font-weight:400; font-style:italic;&quot;&gt;&amp;lt;I&amp;gt; Kurzíva &amp;lt;/I&amp;gt;&lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt; font-style:italic;&quot;&gt;&lt;span style=&quot; font-style:normal;&quot;&gt;(poznámka:  entita &amp;amp;copy; sa používa pre znak copyrightu)&lt;/span&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message encoding="UTF-8">
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;© QGIS 2008&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="unfinished"></translation>
    </message>
</context>
<context>
    <name>QgsCustomProjectionDialog</name>
    <message>
        <source>Delete Projection Definition?</source>
        <translation>Vymazať definíciu zobrazenia?</translation>
    </message>
    <message>
        <source>Deleting a projection definition is not reversable. Do you want to delete it?</source>
        <translation>Vymazanie definície zobrazenia je nevratné. Naozaj si želáte vymazať ju?</translation>
    </message>
    <message>
        <source>Abort</source>
        <translation>Prerušiť</translation>
    </message>
    <message>
        <source>New</source>
        <translation>Nová</translation>
    </message>
    <message>
        <source>QGIS Custom Projection</source>
        <translation>QGIS Vlastné zobrazenia</translation>
    </message>
    <message>
        <source>This proj4 projection definition is not valid. Please correct before pressing save.</source>
        <translation>Táto definícia mapového zobrazenia pre proj4 nie je správna. Pred stlačením tlačidla Uložiť je potrebné ju opraviť.</translation>
    </message>
    <message>
        <source>This proj4 projection definition is not valid. Please give the projection a name before pressing save.</source>
        <translation>Táto definícia mapového zobrazenia pre proj4 nie je správna. Pred stlačením tlačidla Uložiť je potrebné dať zobrazeniu nejaké meno.</translation>
    </message>
    <message>
        <source>This proj4 projection definition is not valid. Please add the parameters before pressing save.</source>
        <translation>Táto definícia mapového zobrazenia pre proj4 nie je správna. Pred stlačením tlačidla Uložiť je potrebné doplniť príslušné parametre zobrazenia.</translation>
    </message>
    <message>
        <source>This proj4 projection definition is not valid. Please add a proj= clause before pressing save.</source>
        <translation>Táto definícia mapového zobrazenia pre proj4 nie je správna. Pred stlačením tlačidla Uložiť je potrebné pridať do definície klauzulu &quot;proj=&quot;.</translation>
    </message>
    <message>
        <source>This proj4 ellipsoid definition is not valid. Please add a ellips= clause before pressing save.</source>
        <translation type="obsolete">Táto definícia elipsoidu pre proj4 nie je správna. Pred stlačením tlačidla Uložiť je potrebné pridať klauzulu &quot;ellips=&quot;.</translation>
    </message>
    <message>
        <source>This proj4 projection definition is not valid.</source>
        <translation>Táto definícia mapového zobrazenia pre proj4 nie je správna.</translation>
    </message>
    <message>
        <source>Northing and Easthing must be in decimal form.</source>
        <translation>Northing (posun na sever) a Easting (posun na východ) musia byť vo forme desatinného čísla.</translation>
    </message>
    <message>
        <source>Internal Error (source projection invalid?)</source>
        <translation>Vnútorná chyba (chybný zdroj mapového zobrazenia?)</translation>
    </message>
</context>
<context>
    <name>QgsCustomProjectionDialogBase</name>
    <message>
        <source>Custom Projection Definition</source>
        <translation type="obsolete">Definícia vlastného zobrazenia</translation>
    </message>
    <message>
        <source>|&lt;</source>
        <translation>|&lt;</translation>
    </message>
    <message>
        <source>&lt;</source>
        <translation>&lt;</translation>
    </message>
    <message>
        <source>1 of 1</source>
        <translation>1 z 1</translation>
    </message>
    <message>
        <source>&gt;</source>
        <translation>&gt;</translation>
    </message>
    <message>
        <source>&gt;|</source>
        <translation>&gt;|</translation>
    </message>
    <message>
        <source>Define</source>
        <translation>Definícia</translation>
    </message>
    <message>
        <source>You can define your own custom projection here. The definition must conform to the proj4 format for specifying a Spatial Reference System.</source>
        <translation type="obsolete">Tu je možné definovať vlastné (mapové) zobrazenie. Definícia musí zodpovedať proj4 formátu špecifikácie priestorového súradnicového systému.</translation>
    </message>
    <message>
        <source>Test</source>
        <translation>Skúška</translation>
    </message>
    <message>
        <source>Calculate</source>
        <translation>Vypočítať</translation>
    </message>
    <message>
        <source>Geographic / WGS84</source>
        <translation>Zemepisný / WGS84</translation>
    </message>
    <message>
        <source>Use the text boxes below to test the projection definition you are creating. Enter a coordinate where both the lat/long and the projected result are known (for example by reading off a map). Then press the calculate button to see if the projection definition you are creating is accurate.</source>
        <translation type="obsolete">Textové polia dole slúžia na otestovanie definície vytvoreného zobrazenia. Vložte súradnice, ktorých hodnoty v danom zobrazení poznáte (napríklad odčítaním z mapy). Následne stlačením tlačitka Vypočítať, možno overiť, či daná definícia zobrazenia je zadaná správne.</translation>
    </message>
    <message>
        <source>Projected Coordinate System</source>
        <translation type="obsolete">Výstupný súradnicový systém</translation>
    </message>
    <message>
        <source>Name</source>
        <translation>Meno</translation>
    </message>
    <message>
        <source>Parameters</source>
        <translation>Parametre</translation>
    </message>
    <message>
        <source>*</source>
        <translation>*</translation>
    </message>
    <message>
        <source>S</source>
        <translation>U</translation>
    </message>
    <message>
        <source>X</source>
        <translation>X</translation>
    </message>
    <message>
        <source>North</source>
        <translation type="unfinished">Sever</translation>
    </message>
    <message>
        <source>East</source>
        <translation type="unfinished">Východ</translation>
    </message>
    <message>
        <source>Custom Coordinate Reference System Definition</source>
        <translation>Definícia vlastného referenčného súradnicového systému</translation>
    </message>
    <message>
        <source>You can define your own custom Coordinate Reference System (CRS) here. The definition must conform to the proj4 format for specifying a CRS.</source>
        <translation>Tu je možné definovať vlastný referenčný súradnicový systém (CRS). Definícia musí zodpovedať proj4 formátu pre špecifikáciu referenčného súradnicového systému.</translation>
    </message>
    <message>
        <source>Use the text boxes below to test the CRS definition you are creating. Enter a coordinate where both the lat/long and the transformed result are known (for example by reading off a map). Then press the calculate button to see if the CRS definition you are creating is accurate.</source>
        <translation>Textové polia dole slúžia na otestovanie definície referenčného súradnicového systému. Zadajte súradnice, ktorých hodnoty v danom súradnicovom systéme poznáte (napríklad odčítaním z mapy). Následne stlačením tlačidla Vypočítať, možno overiť, či je daná definícia súradnicového systému zadaná správne.</translation>
    </message>
    <message>
        <source>Destination CRS        </source>
        <translation>Cieľový ref. súrad. systém</translation>
    </message>
</context>
<context>
    <name>QgsDbSourceSelect</name>
    <message>
        <source>Are you sure you want to remove the </source>
        <translation>Ste si istý, že chcete odstrániť spojenie </translation>
    </message>
    <message>
        <source> connection and all associated settings?</source>
        <translation> a všetky s ním súvisiace nastavenia?</translation>
    </message>
    <message>
        <source>Confirm Delete</source>
        <translation>Potvrdenie mazania</translation>
    </message>
    <message>
        <source>Select Table</source>
        <translation>Vyberte tabuľku</translation>
    </message>
    <message>
        <source>You must select a table in order to add a Layer.</source>
        <translation>Aby bolo možné pridať vrstvu, je treba najskôr vybrať tabuľku.</translation>
    </message>
    <message>
        <source>Password for </source>
        <translation>Heslo pre</translation>
    </message>
    <message>
        <source>Please enter your password:</source>
        <translation>Prosím, vložte vaše heslo:</translation>
    </message>
    <message>
        <source>Connection failed</source>
        <translation>Spojenie zlyhalo</translation>
    </message>
    <message>
        <source>Type</source>
        <translation>Typ</translation>
    </message>
    <message>
        <source>Sql</source>
        <translation>Sql</translation>
    </message>
    <message>
        <source>Connection to %1 on %2 failed. Either the database is down or your settings are incorrect.%3Check your username and password and try again.%4The database said:%5%6</source>
        <translation>Spojenie k databáze %1 na %2 zlyhalo. Databáza je mimo prevádzky, alebo sú nesprávne vaše nastavenia.%3Skontrolujte meno používateľa a heslo a skúste znova.%4Výstup z databázy:%5%6</translation>
    </message>
    <message>
        <source>Wildcard</source>
        <translation type="unfinished">Zástupný znak</translation>
    </message>
    <message>
        <source>RegExp</source>
        <translation type="unfinished">RegVyraz</translation>
    </message>
    <message>
        <source>All</source>
        <translation type="unfinished">Všetko</translation>
    </message>
    <message>
        <source>Schema</source>
        <translation type="unfinished">Schéma</translation>
    </message>
    <message>
        <source>Table</source>
        <translation type="unfinished">Tabuľka</translation>
    </message>
    <message>
        <source>Geometry column</source>
        <translation type="unfinished">Stĺpec s geometriou</translation>
    </message>
    <message>
        <source>Accessible tables could not be determined</source>
        <translation type="unfinished">Nemožno určiť dostupné tabuľky</translation>
    </message>
    <message>
        <source>Database connection was successful, but the accessible tables could not be determined.

The error message from the database was:
%1
</source>
        <translation type="unfinished">Spojenie k databáze bolo úspešné, ale nemožno určiť dostupné tabuľky.

Chybové hlásenie databázy:
%1
</translation>
    </message>
    <message>
        <source>No accessible tables found</source>
        <translation type="unfinished">Nenašli sa žiadne dostupné tabuľky</translation>
    </message>
    <message>
        <source>Database connection was successful, but no accessible tables were found.

Please verify that you have SELECT privilege on a table carrying PostGIS
geometry.</source>
        <translation type="unfinished">Spojenie k databáze bolo úspešné, ale nenašli sa žiadne dostupné tabuľky.

Prosím skontrolujte, či máte práva na SELECT pre tabuľku obsahujúcu PostGIS
geomteriu.</translation>
    </message>
</context>
<context>
    <name>QgsDbSourceSelectBase</name>
    <message>
        <source>Add PostGIS Table(s)</source>
        <translation>Pridať PostGIS tabuľku (tabuľky)</translation>
    </message>
    <message>
        <source>Add</source>
        <translation>Pridať</translation>
    </message>
    <message>
        <source>Help</source>
        <translation>Pomocník</translation>
    </message>
    <message>
        <source>F1</source>
        <translation>F1</translation>
    </message>
    <message>
        <source>Connect</source>
        <translation>Spojiť</translation>
    </message>
    <message>
        <source>New</source>
        <translation>Nové</translation>
    </message>
    <message>
        <source>Edit</source>
        <translation>Upraviť</translation>
    </message>
    <message>
        <source>Delete</source>
        <translation>Vymazať</translation>
    </message>
    <message>
        <source>Close</source>
        <translation>Zatvoriť</translation>
    </message>
    <message>
        <source>PostgreSQL Connections</source>
        <translation>PostgreSQL spojenia</translation>
    </message>
    <message>
        <source>Search:</source>
        <translation type="unfinished">Hľadať:</translation>
    </message>
    <message>
        <source>Search mode:</source>
        <translation type="unfinished">Režim vyhľadávania:</translation>
    </message>
    <message>
        <source>Search in columns:</source>
        <translation type="unfinished">Hľadať v stĺpcoch:</translation>
    </message>
    <message>
        <source>Search options...</source>
        <translation type="unfinished">Voľby hľadania...</translation>
    </message>
</context>
<context>
    <name>QgsDbTableModel</name>
    <message>
        <source>Schema</source>
        <translation>Schéma</translation>
    </message>
    <message>
        <source>Table</source>
        <translation>Tabuľka</translation>
    </message>
    <message>
        <source>Type</source>
        <translation>Typ</translation>
    </message>
    <message>
        <source>Geometry column</source>
        <translation>Stĺpec s geometriou</translation>
    </message>
    <message>
        <source>Sql</source>
        <translation>Sql</translation>
    </message>
    <message>
        <source>Point</source>
        <translation>Bod</translation>
    </message>
    <message>
        <source>Multipoint</source>
        <translation>Multibod</translation>
    </message>
    <message>
        <source>Line</source>
        <translation>Línia</translation>
    </message>
    <message>
        <source>Multiline</source>
        <translation>Multilínia</translation>
    </message>
    <message>
        <source>Polygon</source>
        <translation>Polygón</translation>
    </message>
    <message>
        <source>Multipolygon</source>
        <translation>Multipolygón</translation>
    </message>
</context>
<context>
    <name>QgsDelAttrDialogBase</name>
    <message>
        <source>Delete Attributes</source>
        <translation>Zmazať atribúty</translation>
    </message>
    <message>
        <source>OK</source>
        <translation type="obsolete">OK</translation>
    </message>
    <message>
        <source>Cancel</source>
        <translation type="obsolete">Zrušiť</translation>
    </message>
</context>
<context>
    <name>QgsDelimitedTextPlugin</name>
    <message>
        <source>&amp;Delimited text</source>
        <translation>&amp;Oddelený text</translation>
    </message>
    <message>
        <source>&amp;Add Delimited Text Layer</source>
        <translation>&amp;Pridať vrstvu s oddeleným textom</translation>
    </message>
    <message>
        <source>Add a delimited text file as a map layer. </source>
        <translation> Pridá súbor s oddeleným textom ako mapovú vrstvu.</translation>
    </message>
    <message>
        <source>The file must have a header row containing the field names. </source>
        <translation> Súbor musí obsahovať hlavičkový riadok s názvami jednotlivých polí.</translation>
    </message>
    <message>
        <source>X and Y fields are required and must contain coordinates in decimal units.</source>
        <translation>Požadované sú polia X a Y. Tie musia obsahovať súradnice v desiatkových jednotkách.</translation>
    </message>
    <message>
        <source>DelimitedTextLayer</source>
        <translation>Vrstva z oddeleného textu</translation>
    </message>
</context>
<context>
    <name>QgsDelimitedTextPluginGui</name>
    <message>
        <source>No layer name</source>
        <translation>Žiadne meno vrstvy</translation>
    </message>
    <message>
        <source>Please enter a layer name before adding the layer to the map</source>
        <translation>Vložte prosím meno vrstvy skôr, než ju pridáte do mapy </translation>
    </message>
    <message>
        <source>No delimiter</source>
        <translation>Žiadny oddeľovač</translation>
    </message>
    <message>
        <source>Please specify a delimiter prior to parsing the file</source>
        <translation>Prosím uveďte oddeľovač skôr než sa bude preskúmavať súbor</translation>
    </message>
    <message>
        <source>Choose a delimited text file to open</source>
        <translation>Vyberte súbor s oddeleným textom</translation>
    </message>
    <message>
        <source>Parse</source>
        <translation>Preskúmať</translation>
    </message>
    <message>
        <source>Description</source>
        <translation>Popis</translation>
    </message>
    <message>
        <source>Select a delimited text file containing a header row and one or more rows of x and y coordinates that you would like to use as a point layer and this plugin will do the job for you!</source>
        <translation type="unfinished">Vyberte súbor s oddelným textom obsahujúcim x-ové a y-ové súradnice, ktorý chcete použiť ako vrstvu bodov a tento zásuvný modul sa o zvyšok postará!</translation>
    </message>
    <message>
        <source>Use the layer name box to specify the legend name for the new layer. Use the delimiter box to specify what delimeter is used in your file (e.g. space, comma, tab or a regular expression in Perl style). After choosing a delimiter, press the parse button and select the columns containing the x and y values for the layer.</source>
        <translation type="unfinished">Použite políčko s menom vrstvy na určenie mena novej vrstvy v legende. Políčkom Oddeľovač určíte oddeľovač použitý vo vašom súbore (napr. medzera, čiarka, tabulátor alebo regulárny výraz podobne ako v jazyku Perl). Po vybratí oddeľovača, kliknite na tlačítko Preskúmať a vyberte stĺpce obsahujúce x-ové a a y-ové hodnoty vrstvy.</translation>
    </message>
</context>
<context>
    <name>QgsDelimitedTextPluginGuiBase</name>
    <message>
        <source>Create a Layer from a Delimited Text File</source>
        <translation>Tvorba vrstvy z textového súboru s oddeleným textom (Delimited Text File)</translation>
    </message>
    <message>
        <source>&lt;p align=&quot;right&quot;&gt;X field&lt;/p&gt;</source>
        <translation>&lt;p align=&quot;right&quot;&gt;Pole X&lt;/p&gt;</translation>
    </message>
    <message>
        <source>Name of the field containing x values</source>
        <translation>Meno poľa obsahujúceho x-ové hodnoty</translation>
    </message>
    <message>
        <source>Name of the field containing x values. Choose a field from the list. The list is generated by parsing the header row of the delimited text file.</source>
        <translation>Meno poľa obsahujúceho x-ové hodnoty. Pole vyberte zo zoznamu. Zoznam je vytvorený na základe preskúmania hlavičkového riadku súboru s oddeleným textom.</translation>
    </message>
    <message>
        <source>&lt;p align=&quot;right&quot;&gt;Y field&lt;/p&gt;</source>
        <translation>&lt;p align=&quot;right&quot;&gt;Pole Y&lt;/p&gt;</translation>
    </message>
    <message>
        <source>Name of the field containing y values</source>
        <translation>Meno poľa obsahujúceho y-ové hodnoty</translation>
    </message>
    <message>
        <source>Name of the field containing y values. Choose a field from the list. The list is generated by parsing the header row of the delimited text file.</source>
        <translation>Meno poľa obsahujúceho y-ové hodnoty. Pole vyberte zo zoznamu. Zoznam je generovaný na základe preskúmania hlavičkového riadku súboru s oddeleným textom.</translation>
    </message>
    <message>
        <source>Layer name</source>
        <translation>Meno vrstvy</translation>
    </message>
    <message>
        <source>Name to display in the map legend</source>
        <translation>Meno, ktoré sa objaví v mapovej legende</translation>
    </message>
    <message>
        <source>Name displayed in the map legend</source>
        <translation>Názov zobrazený v okne Legenda</translation>
    </message>
    <message>
        <source>Delimiter</source>
        <translation>Oddeľovač</translation>
    </message>
    <message>
        <source>Delimiter to use when splitting fields in the text file. The delimiter can be more than one character.</source>
        <translation>Oddeľovač používaný na oddelenie jednotlivých polí v textovom súbore. Oddeľovač sa môže skladať aj viac než z jedného znaku.</translation>
    </message>
    <message>
        <source>Delimiter to use when splitting fields in the delimited text file. The delimiter can be 1 or more characters in length.</source>
        <translation>Oddeľovač používaný na oddelenie jednotlivých polí v textovom súbore s oddeleným textom. Oddeľovač sa môže skladať z jedného alebo viac znakov.</translation>
    </message>
    <message>
        <source>Delimited Text Layer</source>
        <translation>Vrstva z oddeleného textu</translation>
    </message>
    <message>
        <source>Delimited text file</source>
        <translation>Súbor s oddeleným textom</translation>
    </message>
    <message>
        <source>Full path to the delimited text file</source>
        <translation>Úplná cesta k textovému súboru s oddeleným textom</translation>
    </message>
    <message>
        <source>Full path to the delimited text file. In order to properly parse the fields in the file, the delimiter must be defined prior to entering the file name. Use the Browse button to the right of this field to choose the input file.</source>
        <translation>Úplná cesta k textovému súboru s oddeleným textom (Delimited Text File). Na to, aby bolo možné preskúmať polia v súbore, musí byť najskôr (prv než meno súboru) správne nastavený oddeľovač textu. Na výber vstupného súboru možno použiť tlačidlo napravo od tohoto poľa.</translation>
    </message>
    <message>
        <source>Browse to find the delimited text file to be processed</source>
        <translation>Prechádzať a nájsť súbor s oddeleným textom, ktorý bude spracovaný</translation>
    </message>
    <message>
        <source>Use this button to browse to the location of the delimited text file. This button will not be enabled until a delimiter has been entered in the &lt;i&gt;Delimiter&lt;/i&gt; box. Once a file is chosen, the X and Y field drop-down boxes will be populated with the fields from the delimited text file.</source>
        <translation>Toto tlačidlo slúži na nájdenie umiestnenia, kde sa nachádza súbor s oddeleným textom. Tlačidlo nebude prístupné kým nebude vyplnené pole &lt;i&gt;Oddeľovač&lt;/i&gt;. Po vybratí súboru budú automaticky (podľa obsahu súboru s oddeleným textom) naplnené aj rozbaľovacie menu s poliami X a Y.</translation>
    </message>
    <message>
        <source>Sample text</source>
        <translation>Vzorka textu</translation>
    </message>
    <message>
        <source>Browse...</source>
        <translation>Prechádzať...</translation>
    </message>
    <message>
        <source>The delimiter is taken as is</source>
        <translation>Presne zadaný  znak oddeľovača</translation>
    </message>
    <message>
        <source>Plain characters</source>
        <translation>Jednoduché znaky</translation>
    </message>
    <message>
        <source>The delimiter is a regular expression</source>
        <translation>Oddeľovačom je regulárny výraz</translation>
    </message>
    <message>
        <source>Regular expression</source>
        <translation>Regulárny výraz</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot;-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="unfinished">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot;-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
</context>
<context>
    <name>QgsDelimitedTextProvider</name>
    <message>
        <source>Note: the following lines were not loaded because Qgis was unable to determine values for the x and y coordinates:
</source>
        <translation>Poznámka: nasledujúce riadky neboli nahraté pretože QGIS nedokázal rozpoznať hodnoty pre súradnice x a y:
</translation>
    </message>
    <message>
        <source>Error</source>
        <translation>Chyba</translation>
    </message>
</context>
<context>
    <name>QgsDetailedItemWidgetBase</name>
    <message>
        <source>Form</source>
        <translation type="unfinished">Vlastnosti mapovej položky</translation>
    </message>
    <message>
        <source>Heading Label</source>
        <translation type="unfinished">Záhlavie</translation>
    </message>
    <message>
        <source>Detail label</source>
        <translation type="unfinished">Detail</translation>
    </message>
</context>
<context>
    <name>QgsDlgPgBufferBase</name>
    <message>
        <source>Buffer features</source>
        <translation>Tvorba okolia (buffera) objektov</translation>
    </message>
    <message>
        <source>Buffer distance in map units:</source>
        <translation>Veľkosť okolia (buffera) v mapových jednotkách:</translation>
    </message>
    <message>
        <source>Table name for the buffered layer:</source>
        <translation>Meno tabuľky pre vrstvu s okolím (bufferom):</translation>
    </message>
    <message>
        <source>Create unique object id</source>
        <translation>Vytvoriť objekt s jedinečným id</translation>
    </message>
    <message>
        <source>public</source>
        <translation>public</translation>
    </message>
    <message>
        <source>Geometry column:</source>
        <translation>Stĺpec s geometriou:</translation>
    </message>
    <message>
        <source>Spatial reference ID:</source>
        <translation>ID priestorového ref. systému:</translation>
    </message>
    <message>
        <source>Unique field to use as feature id:</source>
        <translation>Jedinečné pole použiť ako id objektu:</translation>
    </message>
    <message>
        <source>Schema:</source>
        <translation>Schéma:</translation>
    </message>
    <message>
        <source>Add the buffered layer to the map?</source>
        <translation>Pridať do mapy vrstvu s okolím (bufferom)?</translation>
    </message>
    <message>
        <source>&lt;h2&gt;Buffer the features in layer: &lt;/h2&gt;</source>
        <translation>&lt;h2&gt;Vytvoriť okolie (buffer) objektov vrstvy: &lt;/h2&gt;</translation>
    </message>
    <message>
        <source>Parameters</source>
        <translation>Parametre</translation>
    </message>
</context>
<context>
    <name>QgsEditReservedWordsBase</name>
    <message>
        <source>Edit Reserved Words</source>
        <translation type="obsolete">Upraviť rezervované slovál</translation>
    </message>
    <message>
        <source>Status</source>
        <translation type="obsolete">Stav</translation>
    </message>
    <message>
        <source>Index</source>
        <translation type="obsolete">Index</translation>
    </message>
    <message>
        <source>Reserved Words</source>
        <translation type="obsolete">Rezervované slová</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;/head&gt;&lt;body style=&quot; white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Double click the Column Name column to change the name of the column.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;/head&gt;&lt;body style=&quot; white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Dvojklikom na meno stĺpca je možné zmeniť meno tohoto stĺpca.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>Column Name</source>
        <translation type="obsolete">Meno stĺpca</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;/head&gt;&lt;body style=&quot; white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;This shapefile contains reserved words. These may affect the import into PostgreSQL. Edit the column names so none of the reserved words listed at the right are used (click on a Column Name entry to edit). You may also change any other column name if desired.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;/head&gt;&lt;body style=&quot; white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Tento súbor shape obsahuje rezervované slová. Tieto môžu mať vplyv na import do databázy PostgreSQL. Upravte mená stĺpcov tak, aby nebolo použité žiadne z rezervovaných slov uvedených napravo (upraviť ho možno kliknitím na položku Meno stĺpca). Taktiež možno zmeniť meno aj ľubovoľného ďalšieho stĺpca.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
</context>
<context>
    <name>QgsEncodingFileDialog</name>
    <message>
        <source>Encoding:</source>
        <translation>Kódovanie:</translation>
    </message>
</context>
<context>
    <name>QgsFillStyleWidgetBase</name>
    <message>
        <source>Form1</source>
        <translation type="obsolete">Štýl výplne</translation>
    </message>
    <message>
        <source>Fill Style</source>
        <translation type="obsolete">Štýl výplne</translation>
    </message>
    <message>
        <source>PolyStyleWidget</source>
        <translation type="obsolete">PolyStyleWidget</translation>
    </message>
    <message>
        <source>Colour:</source>
        <translation type="obsolete">Farba:</translation>
    </message>
    <message>
        <source>col</source>
        <translation type="obsolete">col</translation>
    </message>
</context>
<context>
    <name>QgsGPSDeviceDialog</name>
    <message>
        <source>New device %1</source>
        <translation>Nové zariadenie %1</translation>
    </message>
    <message>
        <source>Are you sure?</source>
        <translation>Ste si istý?</translation>
    </message>
    <message>
        <source>Are you sure that you want to delete this device?</source>
        <translation>Ste si istý, že chcete vymazať toto zariadenie?</translation>
    </message>
</context>
<context>
    <name>QgsGPSDeviceDialogBase</name>
    <message>
        <source>GPS Device Editor</source>
        <translation>Editor zariadení GPS</translation>
    </message>
    <message>
        <source>Device name:</source>
        <translation type="obsolete">Názov zariadenia:</translation>
    </message>
    <message>
        <source>This is the name of the device as it will appear in the lists</source>
        <translation>Názov zariadenia tak, ako sa bude objavovať v zoznamoch</translation>
    </message>
    <message>
        <source>Update device</source>
        <translation>Aktualizovať zariadenie</translation>
    </message>
    <message>
        <source>Delete device</source>
        <translation>Vymazať zariadenie</translation>
    </message>
    <message>
        <source>New device</source>
        <translation>Nové zariadenie</translation>
    </message>
    <message>
        <source>Close</source>
        <translation type="obsolete">Zatvoriť</translation>
    </message>
    <message>
        <source>Commands</source>
        <translation>Príkazy</translation>
    </message>
    <message>
        <source>Waypoint download:</source>
        <translation>Stiahnutie orientačných bodov (waypoints):</translation>
    </message>
    <message>
        <source>Waypoint upload:</source>
        <translation>Nahrávanie orientačných bodov (waypoints):</translation>
    </message>
    <message>
        <source>Route download:</source>
        <translation>Stiahnutie cesty (route):</translation>
    </message>
    <message>
        <source>Route upload:</source>
        <translation>Nahrávanie cesty (route):</translation>
    </message>
    <message>
        <source>Track download:</source>
        <translation>Stiahnutie stopy (track):</translation>
    </message>
    <message>
        <source>The command that is used to upload tracks to the device</source>
        <translation>Príkaz, ktorý sa používa na nahratie stôp (tracks) do zariadenia</translation>
    </message>
    <message>
        <source>Track upload:</source>
        <translation>Nahrávanie stôp (track):</translation>
    </message>
    <message>
        <source>The command that is used to download tracks from the device</source>
        <translation>Príkaz, ktorý sa používa na stiahnutie stôp (tracks) z tohoto zariadenia</translation>
    </message>
    <message>
        <source>The command that is used to upload routes to the device</source>
        <translation>Príkaz, ktorý sa používa na nahratie ciest (rosutes) do tohoto zariadenia</translation>
    </message>
    <message>
        <source>The command that is used to download routes from the device</source>
        <translation>Príkaz, ktorý sa používa na stiahnutie ciest (routes) z tohoto zariadenia</translation>
    </message>
    <message>
        <source>The command that is used to upload waypoints to the device</source>
        <translation>Príkaz, ktorý sa používa na nahratie orientačných bodov (waypoints) do tohoto zariadenia</translation>
    </message>
    <message>
        <source>The command that is used to download waypoints from the device</source>
        <translation>Príkaz, ktorý sa používa na stiahnutie orientačných bodov (waypoints) z tohoto zariadenia</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;/head&gt;&lt;body style=&quot; white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;&lt;p style=&quot; margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;In the download and upload commands there can be special words that will be replaced by QGIS when the commands are used. These words are:&lt;span style=&quot; font-style:italic;&quot;&gt;%babel&lt;/span&gt; - the path to GPSBabel&lt;br /&gt;&lt;span style=&quot; font-style:italic;&quot;&gt;%in&lt;/span&gt; - the GPX filename when uploading or the port when downloading&lt;br /&gt;&lt;span style=&quot; font-style:italic;&quot;&gt;%out&lt;/span&gt; - the port when uploading or the GPX filename when downloading&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;/head&gt;&lt;body style=&quot; white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;&lt;p style=&quot; margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;p&gt;V príkazoch na stiahnutie a nahratie, môžu byť použité špeciálne reťazce, ktoré budú pri použití QGISom nahradené. Sú to:&lt;span style=&quot; font-style:italic;&quot;&gt;%babel&lt;/span&gt; - cesta k GPSBabel&lt;br /&gt;&lt;span style=&quot; font-style:italic;&quot;&gt;%in&lt;/span&gt; - meno súboru GPX pri nahrávaní (zo zariadenia GPS), alebo port pri nahrávaní (do zariadenia GPS)&lt;br /&gt;&lt;span style=&quot; font-style:italic;&quot;&gt;%out&lt;/span&gt; - port pri nahrávaní (do zariadenia GPS), alebo meno súboru GPX pri nahrávaní (zo zariadenia GPS)&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>Device name</source>
        <translation>Názov zariadenia</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;In the download and upload commands there can be special words that will be replaced by QGIS when the commands are used. These words are:&lt;span style=&quot; font-style:italic;&quot;&gt;%babel&lt;/span&gt; - the path to GPSBabel&lt;br /&gt;&lt;span style=&quot; font-style:italic;&quot;&gt;%in&lt;/span&gt; - the GPX filename when uploading or the port when downloading&lt;br /&gt;&lt;span style=&quot; font-style:italic;&quot;&gt;%out&lt;/span&gt; - the port when uploading or the GPX filename when downloading&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="unfinished">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;V príkazoch na stiahnutie a nahratie, môžu byť použité špeciálne reťazce, ktoré budú pri použití QGISom nahradené. Sú to:&lt;span style=&quot; font-style:italic;&quot;&gt;%babel&lt;/span&gt; - cesta k GPSBabel&lt;br /&gt;&lt;span style=&quot; font-style:italic;&quot;&gt;%in&lt;/span&gt;  - meno súboru GPX pri nahrávaní (zo zariadenia GPS), alebo port pri nahrávaní (do zariadenia GPS)&lt;br /&gt;&lt;span style=&quot; font-style:italic;&quot;&gt;%out&lt;/span&gt; port pri nahrávaní (do zariadenia GPS), alebo meno súboru GPX pri nahrávaní (zo zariadenia GPS)&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
</context>
<context>
    <name>QgsGPSPlugin</name>
    <message>
        <source>&amp;Gps</source>
        <translation>&amp;GPS</translation>
    </message>
    <message>
        <source>&amp;Gps Tools</source>
        <translation>&amp;Nástroje GPS</translation>
    </message>
    <message>
        <source>&amp;Create new GPX layer</source>
        <translation>Vytvoriť &amp;novú vrstvu GPX</translation>
    </message>
    <message>
        <source>Creates a new GPX layer and displays it on the map canvas</source>
        <translation>Vytvorí novú vrstvu GPX a zobrazí ju v mapovom pohľade</translation>
    </message>
    <message>
        <source>Save new GPX file as...</source>
        <translation>Uložiť novú vrstvu GPX ako...</translation>
    </message>
    <message>
        <source>GPS eXchange file (*.gpx)</source>
        <translation>Výmennný súbor GPS (*.gpx)</translation>
    </message>
    <message>
        <source>Could not create file</source>
        <translation>Nemožno vytvoriť súbor</translation>
    </message>
    <message>
        <source>Unable to create a GPX file with the given name. </source>
        <translation> Nemožno vytvoriť súbor GPX s daným menom.</translation>
    </message>
    <message>
        <source>Try again with another name or in another </source>
        <translation> Zmeňte meno súboru alebo adresár a skúste</translation>
    </message>
    <message>
        <source>directory.</source>
        <translation>znova.</translation>
    </message>
    <message>
        <source>GPX Loader</source>
        <translation>Nahrávač GPX</translation>
    </message>
    <message>
        <source>Unable to read the selected file.
</source>
        <translation>Vybraný súbor nemožno prečítať.
</translation>
    </message>
    <message>
        <source>Please reselect a valid file.</source>
        <translation>Prosím vyberte platný súbor.</translation>
    </message>
    <message>
        <source>Could not start process</source>
        <translation>Nemožno spustiť proces</translation>
    </message>
    <message>
        <source>Could not start GPSBabel!</source>
        <translation>Nemožno spustiť GPSBabel!</translation>
    </message>
    <message>
        <source>Importing data...</source>
        <translation>Importujú sa údaje...</translation>
    </message>
    <message>
        <source>Cancel</source>
        <translation>Zrušiť</translation>
    </message>
    <message>
        <source>Could not import data from %1!

</source>
        <translation>Nemožno importovať údaje z %1!
</translation>
    </message>
    <message>
        <source>Error importing data</source>
        <translation>Chyba pri importe údajov</translation>
    </message>
    <message>
        <source>Not supported</source>
        <translation>Nepodporované</translation>
    </message>
    <message>
        <source>This device does not support downloading </source>
        <translation> Toto zariadenie nepodporuje stiahnutie (download) údajov</translation>
    </message>
    <message>
        <source>of </source>
        <translation> z</translation>
    </message>
    <message>
        <source>Downloading data...</source>
        <translation>Sťahujú sa údaje...</translation>
    </message>
    <message>
        <source>Could not download data from GPS!

</source>
        <translation>Nemožno stiahnuť údaje z GPS!
</translation>
    </message>
    <message>
        <source>Error downloading data</source>
        <translation>Chyba pri sťahovaní údajov</translation>
    </message>
    <message>
        <source>This device does not support uploading of </source>
        <translation> Toto zariadenie nepodporuje nahrávanie (upload) z </translation>
    </message>
    <message>
        <source>Uploading data...</source>
        <translation>Nahrávajú sa údaje...</translation>
    </message>
    <message>
        <source>Error while uploading data to GPS!

</source>
        <translation>Chyba pri nahrávaní údajov do GPS!
</translation>
    </message>
    <message>
        <source>Error uploading data</source>
        <translation>Chyba pri nahrávaní údajov</translation>
    </message>
    <message>
        <source>Could not convert data from %1!

</source>
        <translation>Nemožno previesť údaje z %1!

</translation>
    </message>
    <message>
        <source>Error converting data</source>
        <translation>Chyba pri prevode údajov</translation>
    </message>
</context>
<context>
    <name>QgsGPSPluginGui</name>
    <message>
        <source>Choose a filename to save under</source>
        <translation type="obsolete">Vyberte meno súboru do ktorého sa bude ukladať</translation>
    </message>
    <message>
        <source>GPS eXchange format (*.gpx)</source>
        <translation>Výmennný súbor GPS (*.gpx)</translation>
    </message>
    <message>
        <source>Select GPX file</source>
        <translation>Vyberte súbor GPX</translation>
    </message>
    <message>
        <source>Select file and format to import</source>
        <translation>Vyberte súbor a formát na import</translation>
    </message>
    <message>
        <source>Waypoints</source>
        <translation>Orientačné body (waypoints)</translation>
    </message>
    <message>
        <source>Routes</source>
        <translation>Cesty (routes)</translation>
    </message>
    <message>
        <source>Tracks</source>
        <translation>Stopy (tracks)</translation>
    </message>
    <message>
        <source>QGIS can perform conversions of GPX files, by using GPSBabel (%1) to perform the conversions.</source>
        <translation type="unfinished">QGIS dokáže robiť prevody GPX súborov, na ktoré využíva GPSBabel (%1).</translation>
    </message>
    <message>
        <source>This requires that you have GPSBabel installed where QGIS can find it.</source>
        <translation type="unfinished">To vyžaduje aby bol GPSBabel nainštalovaný tam kde ho môže QGIS nájsť.</translation>
    </message>
    <message>
        <source>Select a GPX input file name, the type of conversion you want to perform, a GPX filename that you want to save the converted file as, and a name for the new layer created from the result.</source>
        <translation type="obsolete">Vyberte meno vstupného súboru GPX, typ prevodu ktorý chcete vykonať a meno súboru GPX do ktorého bude uložený výsledný súbor a meno novej výslednej vrstvy.</translation>
    </message>
    <message>
        <source>GPX is the %1, which is used to store information about waypoints, routes, and tracks.</source>
        <translation type="unfinished">GPX je %1 a používa sa na ukladanie informácií o orientačných bodoch, cestách a stopách.</translation>
    </message>
    <message>
        <source>GPS eXchange file format</source>
        <translation type="unfinished">výmenným formátom údajov z GPS (GPS eXchange)</translation>
    </message>
    <message>
        <source>Select a GPX file and then select the feature types that you want to load.</source>
        <translation type="unfinished">Vyberte GPX súbor a následne typy objektu, ktoré chcete nahrať.</translation>
    </message>
    <message>
        <source>This tool will help you download data from a GPS device.</source>
        <translation type="unfinished">Tento nástroj vám pomôže stiahnuť údaje zo zariadenia GPS.</translation>
    </message>
    <message>
        <source>Choose your GPS device, the port it is connected to, the feature type you want to download, a name for your new layer, and the GPX file where you want to store the data.</source>
        <translation type="unfinished">Vyberte zariadenie GPS, port ku ktorému je pripojené, typ objektov, ktorý sa má sťahovať, meno novej vrstvy, a GPX súboru, kam sa majú uložiť údaje.</translation>
    </message>
    <message>
        <source>If your device isn&apos;t listed, or if you want to change some settings, you can also edit the devices.</source>
        <translation type="unfinished">Ak vaše zariadenie ne je v zozname, alebo ak je potrebné zmeniť niektoré nastavenia, je možné tiež upravovať vlastnosti jednotlivých zariadení.</translation>
    </message>
    <message>
        <source>This tool uses the program GPSBabel (%1) to transfer the data.</source>
        <translation type="unfinished">Tento nástroj využíva na prenos údajov program GPSBabel (%1).</translation>
    </message>
    <message>
        <source>This tool will help you upload data from a GPX layer to a GPS device.</source>
        <translation type="unfinished">Tento nástroj vám pomôže nahrať údaje do zariadenia GPS.</translation>
    </message>
    <message>
        <source>Choose the layer you want to upload, the device you want to upload it to, and the port your device is connected to.</source>
        <translation type="unfinished">Vyberte vrstvu, ktorú chcete nahrať, zariadenie do ktorého ju chcete nahrať, a port, ku ktorému je zariadenie pripojené.</translation>
    </message>
    <message>
        <source>QGIS can only load GPX files by itself, but many other formats can be converted to GPX using GPSBabel (%1).</source>
        <translation type="unfinished">QGIS môže nahrať GPX súbory ako také, ale mnoho ďalších formátov je možné skonvertovať do GPX s použitím nástroja GPSBabel (%1). </translation>
    </message>
    <message>
        <source>Select a GPS file format and the file that you want to import, the feature type that you want to use, a GPX filename that you want to save the converted file as, and a name for the new layer.</source>
        <translation type="obsolete">Vyberte formát súboru s údajmi z GPS a súbor ktorý bude importovaný, typ objektu, ktorý chcete použiť, meno súboru GPX, do ktorého bude uložený výsledný súbor a meno novej vrstvy.</translation>
    </message>
    <message>
        <source>All file formats can not store waypoints, routes, and tracks, so some feature types may be disabled for some file formats.</source>
        <translation type="unfinished">Vo všetkých súborových formátoch nie je možné mať uložené orientačné body, cesty, či stopy, takže niektoré typy objektov môžu byť pre určité formáty súborov zakázané.</translation>
    </message>
    <message>
        <source>Choose a file name to save under</source>
        <translation type="unfinished">Vyberte meno súboru do ktorého sa bude ukladať</translation>
    </message>
    <message>
        <source>Select a GPS file format and the file that you want to import, the feature type that you want to use, a GPX file name that you want to save the converted file as, and a name for the new layer.</source>
        <translation type="unfinished">Vyberte formát súboru s údajmi z GPS a súbor ktorý bude importovaný, typ objektu, ktorý chcete použiť, meno súboru GPX, do ktorého bude uložený výsledný súbor a meno novej vrstvy.</translation>
    </message>
    <message>
        <source>Select a GPX input file name, the type of conversion you want to perform, a GPX file name that you want to save the converted file as, and a name for the new layer created from the result.</source>
        <translation type="unfinished">Vyberte meno vstupného súboru GPX, typ prevodu ktorý chcete vykonať a meno súboru GPX do ktorého bude uložený výsledný súbor a meno novej výslednej vrstvy.</translation>
    </message>
</context>
<context>
    <name>QgsGPSPluginGuiBase</name>
    <message>
        <source>GPS Tools</source>
        <translation>Nástroje na prácu s GPS</translation>
    </message>
    <message>
        <source>Load GPX file</source>
        <translation>Nahrať GPX súbor</translation>
    </message>
    <message>
        <source>File:</source>
        <translation>Súbor:</translation>
    </message>
    <message>
        <source>Feature types:</source>
        <translation>Typy objektov:</translation>
    </message>
    <message>
        <source>Waypoints</source>
        <translation>Orientačné body (waypoints)</translation>
    </message>
    <message>
        <source>Routes</source>
        <translation>Cesty (routes)</translation>
    </message>
    <message>
        <source>Tracks</source>
        <translation>Stopy (tracks)</translation>
    </message>
    <message>
        <source>Import other file</source>
        <translation>Importovať iný súbor</translation>
    </message>
    <message>
        <source>File to import:</source>
        <translation>Súbor na import:</translation>
    </message>
    <message>
        <source>Feature type:</source>
        <translation>Typ objektu:</translation>
    </message>
    <message>
        <source>GPX output file:</source>
        <translation>Výstupný súbor GPX:</translation>
    </message>
    <message>
        <source>Layer name:</source>
        <translation>Názov vrstvy:</translation>
    </message>
    <message>
        <source>Download from GPS</source>
        <translation>Stiahnuť z GPS</translation>
    </message>
    <message>
        <source>Edit devices</source>
        <translation>Upraviť zariadenia</translation>
    </message>
    <message>
        <source>GPS device:</source>
        <translation>GPS zariadenie:</translation>
    </message>
    <message>
        <source>Output file:</source>
        <translation>Výstupný súbor:</translation>
    </message>
    <message>
        <source>Port:</source>
        <translation>Port:</translation>
    </message>
    <message>
        <source>Upload to GPS</source>
        <translation>Nahrať do GPS</translation>
    </message>
    <message>
        <source>Data layer:</source>
        <translation>Vrstva údajov:</translation>
    </message>
    <message>
        <source>Browse...</source>
        <translation>Prechádzať...</translation>
    </message>
    <message>
        <source>Save As...</source>
        <translation>Uložiť ako...</translation>
    </message>
    <message>
        <source>(Note: Selecting correct file type in browser dialog important!)</source>
        <translation type="unfinished">(Poznámka: Je dôležité vybrať v dialógovom okne výberu súboru správny typ súboru!)</translation>
    </message>
    <message>
        <source>GPX Conversions</source>
        <translation type="unfinished">Prevody GPX</translation>
    </message>
    <message>
        <source>Conversion:</source>
        <translation type="unfinished">Prevod:</translation>
    </message>
    <message>
        <source>GPX input file:</source>
        <translation>Vstupný súbor GPX:</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;
&lt;p style=&quot;-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;
&lt;p style=&quot;-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot;-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="unfinished">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot;-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>Edit devices...</source>
        <translation>Upraviť zariadenia...</translation>
    </message>
    <message>
        <source>Refresh</source>
        <translation>Obnoviť</translation>
    </message>
</context>
<context>
    <name>QgsGPXProvider</name>
    <message>
        <source>Bad URI - you need to specify the feature type.</source>
        <translation>Chybné URI - je potrebné manuálne určiť typ objektu.</translation>
    </message>
    <message>
        <source>GPS eXchange file</source>
        <translation>Súbor GPS eXchange</translation>
    </message>
    <message>
        <source>Digitized in QGIS</source>
        <translation>Digitalizované v QGISe</translation>
    </message>
</context>
<context>
    <name>QgsGenericProjectionSelector</name>
    <message>
        <source>Define this layer&apos;s projection:</source>
        <translation>Určiť mapové zobrazenie tejto vrstvy:</translation>
    </message>
    <message>
        <source>This layer appears to have no projection specification.</source>
        <translation>Táto vrstva zrejme nemá uvedenú žiadnu informáciu o použitom mapovom zobrazení.</translation>
    </message>
    <message>
        <source>By default, this layer will now have its projection set to that of the project, but you may override this by selecting a different projection below.</source>
        <translation>Predvolene bude mať teraz nastavené mapové zobrazenie zhodné s nastavením projektu, ale je možné ho zmeniť vybratím iného zobraznia z nižšie uvedenej ponuky.</translation>
    </message>
</context>
<context>
    <name>QgsGenericProjectionSelectorBase</name>
    <message>
        <source>Projection Selector</source>
        <translation>Výber zobrazenia</translation>
    </message>
</context>
<context>
    <name>QgsGeomTypeDialog</name>
    <message>
        <source>Real</source>
        <translation type="unfinished">Reálne číslo</translation>
    </message>
    <message>
        <source>Integer</source>
        <translation type="unfinished">Celé číslo</translation>
    </message>
    <message>
        <source>String</source>
        <translation type="unfinished">Reťazec</translation>
    </message>
</context>
<context>
    <name>QgsGeomTypeDialogBase</name>
    <message>
        <source>Type</source>
        <translation>Typ</translation>
    </message>
    <message>
        <source>Point</source>
        <translation>Bod</translation>
    </message>
    <message>
        <source>Line</source>
        <translation>Línia</translation>
    </message>
    <message>
        <source>Polygon</source>
        <translation>Polygón</translation>
    </message>
    <message>
        <source>New Vector Layer</source>
        <translation>Nová vektorová vrstva</translation>
    </message>
    <message>
        <source>Add</source>
        <translation type="obsolete">Pridať</translation>
    </message>
    <message>
        <source>Remove</source>
        <translation type="obsolete">Odobrať</translation>
    </message>
    <message>
        <source>File format</source>
        <translation>Formát súboru</translation>
    </message>
    <message>
        <source>Attributes</source>
        <translation>Atribúty</translation>
    </message>
    <message>
        <source>Name</source>
        <translation type="unfinished">Meno</translation>
    </message>
    <message>
        <source>Remove selected row</source>
        <translation type="obsolete">Odstrániť vybratý riadok</translation>
    </message>
    <message>
        <source>...</source>
        <translation type="unfinished">...</translation>
    </message>
    <message>
        <source>Add values manually</source>
        <translation type="obsolete">Pridať hodnoty ručne</translation>
    </message>
    <message>
        <source>Delete selected attribute</source>
        <translation type="unfinished">Zmazať vybraný atribút</translation>
    </message>
    <message>
        <source>Add attribute</source>
        <translation>Pridať atribút</translation>
    </message>
</context>
<context>
    <name>QgsGeorefDescriptionDialogBase</name>
    <message>
        <source>Description georeferencer</source>
        <translation>Popis Georeferencera</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:12pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;span style=&quot; font-size:11pt; font-weight:600;&quot;&gt;Description&lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:9pt;&quot;&gt;This plugin can generate world files for rasters. You select points on the raster and give their world coordinates, and the plugin will compute the world file parameters. The more coordinates you can provide the better the result will be.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="unfinished">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:12pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;span style=&quot; font-size:11pt; font-weight:600;&quot;&gt;Popis&lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:9pt;&quot;&gt;Tento zásuvný dokáže vygenerovať world súbory pre rastre. Stačí vybrať body na rastri a určiť ich skutočné world súradnice, a modul vypočítal parametre world súboru. Viac súradníc umožňuje dosiahnuť lepšie výsledky.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
</context>
<context>
    <name>QgsGeorefPlugin</name>
    <message>
        <source>&amp;Georeferencer</source>
        <translation>&amp;Georeferencer</translation>
    </message>
</context>
<context>
    <name>QgsGeorefPluginGui</name>
    <message>
        <source>Choose a raster file</source>
        <translation type="unfinished">Vyberte rastrový súbor</translation>
    </message>
    <message>
        <source>Raster files (*.*)</source>
        <translation>Rastrové súbory (*.*)</translation>
    </message>
    <message>
        <source>Error</source>
        <translation>Chyba</translation>
    </message>
    <message>
        <source>The selected file is not a valid raster file.</source>
        <translation type="unfinished">Vybraný súbor nie je platný rastrový súbor.</translation>
    </message>
    <message>
        <source>World file exists</source>
        <translation>Súbor world existuje</translation>
    </message>
    <message>
        <source>&lt;p&gt;The selected file already seems to have a </source>
        <translation>&lt;p&gt;K vybranému súboru už zrejme existuje</translation>
    </message>
    <message>
        <source>world file! Do you want to replace it with the </source>
        <translation>world súbor! Želáte si nahradiť ho </translation>
    </message>
    <message>
        <source>new world file?&lt;/p&gt;</source>
        <translation>s novým world súborom?&lt;/p&gt;</translation>
    </message>
</context>
<context>
    <name>QgsGeorefPluginGuiBase</name>
    <message>
        <source>Georeferencer</source>
        <translation>Georerenčný modul</translation>
    </message>
    <message>
        <source>Close</source>
        <translation>Zatvoriť</translation>
    </message>
    <message>
        <source>...</source>
        <translation type="unfinished">...</translation>
    </message>
    <message>
        <source>Raster file:</source>
        <translation>Rastrový súbor:</translation>
    </message>
    <message>
        <source>Arrange plugin windows</source>
        <translation type="unfinished">Usporiadať okná zásuvného modulu</translation>
    </message>
    <message>
        <source>Description...</source>
        <translation type="unfinished">Popis...</translation>
    </message>
</context>
<context>
    <name>QgsGeorefWarpOptionsDialog</name>
    <message>
        <source>unstable</source>
        <translation type="obsolete">nestabilná</translation>
    </message>
</context>
<context>
    <name>QgsGeorefWarpOptionsDialogBase</name>
    <message>
        <source>Warp options</source>
        <translation type="unfinished">Voľby zakrivenia</translation>
    </message>
    <message>
        <source>Resampling method:</source>
        <translation>Prevzorkovacia metóda:</translation>
    </message>
    <message>
        <source>Nearest neighbour</source>
        <translation>Metóda najbližšieho suseda</translation>
    </message>
    <message>
        <source>Linear</source>
        <translation>Lineárna</translation>
    </message>
    <message>
        <source>Cubic</source>
        <translation>Kubická</translation>
    </message>
    <message>
        <source>OK</source>
        <translation>OK</translation>
    </message>
    <message>
        <source>Use 0 for transparency when needed</source>
        <translation>Ak bude treba, použiť 0 pre priehľadnosť</translation>
    </message>
    <message>
        <source>Compression:</source>
        <translation type="unfinished">Druh kompresie:</translation>
    </message>
</context>
<context>
    <name>QgsGraduatedSymbolDialog</name>
    <message>
        <source>Empty</source>
        <translation>Prázdny</translation>
    </message>
    <message>
        <source>Equal Interval</source>
        <translation>Rovnaký interval</translation>
    </message>
    <message>
        <source>Quantiles</source>
        <translation>Kvantily</translation>
    </message>
</context>
<context>
    <name>QgsGraduatedSymbolDialogBase</name>
    <message>
        <source>graduated Symbol</source>
        <translation>stupňovaný symbol</translation>
    </message>
    <message>
        <source>Delete class</source>
        <translation>Zmazať triedu</translation>
    </message>
    <message>
        <source>Classify</source>
        <translation type="unfinished">Triediť</translation>
    </message>
    <message>
        <source>Classification field</source>
        <translation type="unfinished">Pole určujúce triedu</translation>
    </message>
    <message>
        <source>Mode</source>
        <translation type="unfinished">Mód</translation>
    </message>
    <message>
        <source>Number of classes</source>
        <translation type="unfinished">Počet tried</translation>
    </message>
</context>
<context>
    <name>QgsGrassAttributes</name>
    <message>
        <source>Warning</source>
        <translation>Upozornenie</translation>
    </message>
    <message>
        <source>Column</source>
        <translation>Stĺpec</translation>
    </message>
    <message>
        <source>Value</source>
        <translation>Hodnota</translation>
    </message>
    <message>
        <source>Type</source>
        <translation>Typ</translation>
    </message>
    <message>
        <source>ERROR</source>
        <translation>CHYBA</translation>
    </message>
    <message>
        <source>OK</source>
        <translation>OK</translation>
    </message>
    <message>
        <source>Layer</source>
        <translation>Vrstva</translation>
    </message>
</context>
<context>
    <name>QgsGrassAttributesBase</name>
    <message>
        <source>GRASS Attributes</source>
        <translation>GRASS - atribúty</translation>
    </message>
    <message>
        <source>Tab 1</source>
        <translation type="unfinished">Tab 1</translation>
    </message>
    <message>
        <source>result</source>
        <translation>výsledok</translation>
    </message>
    <message>
        <source>Update</source>
        <translation>Aktualizovať</translation>
    </message>
    <message>
        <source>Update database record</source>
        <translation>Aktualizovať záznam v databáze</translation>
    </message>
    <message>
        <source>New</source>
        <translation>Nový</translation>
    </message>
    <message>
        <source>Add new category using settings in GRASS Edit toolbox</source>
        <translation>Pridať novú kategóriu s použitím nastavení v nástrojovom paneli GRASS - Úpravy</translation>
    </message>
    <message>
        <source>Delete</source>
        <translation>Vymazať</translation>
    </message>
    <message>
        <source>Delete selected category</source>
        <translation>Vymazať vybranú kategóriu</translation>
    </message>
</context>
<context>
    <name>QgsGrassBrowser</name>
    <message>
        <source>Tools</source>
        <translation>Nástroje</translation>
    </message>
    <message>
        <source>Add selected map to canvas</source>
        <translation>Pridať vybranú mapu na plátno</translation>
    </message>
    <message>
        <source>Delete selected map</source>
        <translation>Vymazať vybranú mapu</translation>
    </message>
    <message>
        <source>Refresh</source>
        <translation>Obnoviť</translation>
    </message>
    <message>
        <source>Copy selected map</source>
        <translation>Kopírovať vybranú mapu</translation>
    </message>
    <message>
        <source>Rename selected map</source>
        <translation>Premenovať vybranú mapu</translation>
    </message>
    <message>
        <source>Set current region to selected map</source>
        <translation>Nastaviť aktuálny región na rozsah vybranej mapy</translation>
    </message>
    <message>
        <source>Warning</source>
        <translation>Upozornenie</translation>
    </message>
    <message>
        <source>Cannot copy map </source>
        <translation> Nemožno skopírovať mapu</translation>
    </message>
    <message>
        <source>&lt;br&gt;command: </source>
        <translation> &lt;br&gt;príkaz:</translation>
    </message>
    <message>
        <source>Cannot rename map </source>
        <translation> Nemožno premenovať mapu</translation>
    </message>
    <message>
        <source>Delete map &lt;b&gt;</source>
        <translation>Odstrániť mapu &lt;b&gt;</translation>
    </message>
    <message>
        <source>Cannot delete map </source>
        <translation> Nemožno odstrániť mapu</translation>
    </message>
    <message>
        <source>Cannot write new region</source>
        <translation>Nemožno zapísať nový región</translation>
    </message>
    <message>
        <source>New name</source>
        <translation>Nové meno</translation>
    </message>
</context>
<context>
    <name>QgsGrassEdit</name>
    <message>
        <source>New point</source>
        <translation>Nový bod</translation>
    </message>
    <message>
        <source>New centroid</source>
        <translation>Nový centroid</translation>
    </message>
    <message>
        <source>Delete vertex</source>
        <translation>Vymazať uzol</translation>
    </message>
    <message>
        <source>Left: </source>
        <translation> Ľavé tlačidlo: </translation>
    </message>
    <message>
        <source>Middle: </source>
        <translation> Stredné tlačidlo: </translation>
    </message>
    <message>
        <source>Edit tools</source>
        <translation>Nástroje na úpravy</translation>
    </message>
    <message>
        <source>New line</source>
        <translation>Nová línia</translation>
    </message>
    <message>
        <source>New boundary</source>
        <translation>Nové ohraničenie</translation>
    </message>
    <message>
        <source>Move vertex</source>
        <translation>Presunúť uzol</translation>
    </message>
    <message>
        <source>Add vertex</source>
        <translation>Pridať uzol</translation>
    </message>
    <message>
        <source>Split line</source>
        <translation>Rozdeliť líniu</translation>
    </message>
    <message>
        <source>Edit attributes</source>
        <translation>Upraviť atribúty</translation>
    </message>
    <message>
        <source>Close</source>
        <translation>Zatvoriť</translation>
    </message>
    <message>
        <source>Move element</source>
        <translation>Posun prvku</translation>
    </message>
    <message>
        <source>Delete element</source>
        <translation>Vymazať prvok</translation>
    </message>
    <message>
        <source>Warning</source>
        <translation>Upozornenie</translation>
    </message>
    <message>
        <source>You are not owner of the mapset, cannot open the vector for editing.</source>
        <translation>Nie ste vlastníkom tohto súboru máp (mapsetu), vektorovú vrstvu nemožno otvoriť pre úpravy.</translation>
    </message>
    <message>
        <source>Cannot open vector for update.</source>
        <translation type="unfinished">Nemožno otvoriť vektorovú mapu pre update.</translation>
    </message>
    <message>
        <source>Info</source>
        <translation type="unfinished">Info</translation>
    </message>
    <message>
        <source>The table was created</source>
        <translation type="unfinished">Tabuľka bola vytvorená</translation>
    </message>
    <message>
        <source>Tool not yet implemented.</source>
        <translation>Nástroj nie je zatiaľ implementovaný.</translation>
    </message>
    <message>
        <source>Cannot check orphan record: </source>
        <translation>Nemožno skontrolovať osirelý záznam: </translation>
    </message>
    <message>
        <source>Orphan record was left in attribute table. &lt;br&gt;Delete the record?</source>
        <translation>V tabuľke atribútov ostal osirelý záznam. &lt;br&gt;Zmazať tento záznam?</translation>
    </message>
    <message>
        <source>Cannot delete orphan record: </source>
        <translation>Nemožno zmazať osirelý záznam: </translation>
    </message>
    <message>
        <source>Cannot describe table for field </source>
        <translation type="unfinished">Nemožno opísať tabuľku pre pole </translation>
    </message>
    <message>
        <source>Background</source>
        <translation>Pozadie</translation>
    </message>
    <message>
        <source>Highlight</source>
        <translation>Zvýraznenie</translation>
    </message>
    <message>
        <source>Dynamic</source>
        <translation type="unfinished">Dynamický</translation>
    </message>
    <message>
        <source>Point</source>
        <translation>Bod</translation>
    </message>
    <message>
        <source>Line</source>
        <translation>Línia</translation>
    </message>
    <message>
        <source>Boundary (no area)</source>
        <translation>Ohraničenie (nie oblasť)</translation>
    </message>
    <message>
        <source>Boundary (1 area)</source>
        <translation>Ohraničenie (1 oblasť)</translation>
    </message>
    <message>
        <source>Boundary (2 areas)</source>
        <translation>Ohraničenie (2 oblasti)</translation>
    </message>
    <message>
        <source>Centroid (in area)</source>
        <translation>Centroid (v oblasti)</translation>
    </message>
    <message>
        <source>Centroid (outside area)</source>
        <translation>Centroid (mimo oblasti)</translation>
    </message>
    <message>
        <source>Centroid (duplicate in area)</source>
        <translation type="unfinished">Centroid (dvojník v oblasti)</translation>
    </message>
    <message>
        <source>Node (1 line)</source>
        <translation type="unfinished">Uzol (1 línia)</translation>
    </message>
    <message>
        <source>Node (2 lines)</source>
        <translation type="unfinished">Uzol (2 línie)</translation>
    </message>
    <message>
        <source>Disp</source>
        <comment>

Column title</comment>
        <translation type="obsolete">Zobraz</translation>
    </message>
    <message>
        <source>Color</source>
        <comment>Column title</comment>
        <translation type="obsolete">Farba</translation>
    </message>
    <message>
        <source>Type</source>
        <comment>Column title</comment>
        <translation type="obsolete">Typ</translation>
    </message>
    <message>
        <source>Index</source>
        <comment>Column title</comment>
        <translation type="obsolete">Index</translation>
    </message>
    <message>
        <source>Column</source>
        <translation type="obsolete">Stĺpec</translation>
    </message>
    <message>
        <source>Type</source>
        <translation type="obsolete">Typ</translation>
    </message>
    <message>
        <source>Length</source>
        <translation type="obsolete">Dĺžka</translation>
    </message>
    <message>
        <source>Next not used</source>
        <translation>Najbližšie nepoužité</translation>
    </message>
    <message>
        <source>Manual entry</source>
        <translation>Ručné vkladanie</translation>
    </message>
    <message>
        <source>No category</source>
        <translation>Bez kategórie</translation>
    </message>
    <message>
        <source>Right: </source>
        <translation type="unfinished">Pravé: </translation>
    </message>
    <message>
        <source>Disp</source>
        <comment>Column title</comment>
        <translation type="obsolete">Zobraz</translation>
    </message>
</context>
<context>
    <name>QgsGrassEditBase</name>
    <message>
        <source>GRASS Edit</source>
        <translation>GRASS - úpravy</translation>
    </message>
    <message>
        <source>Category</source>
        <translation>Kategória</translation>
    </message>
    <message>
        <source>Mode</source>
        <translation>Mód</translation>
    </message>
    <message>
        <source>Settings</source>
        <translation>Nastavenia</translation>
    </message>
    <message>
        <source>Snapping in screen pixels</source>
        <translation>Zameriavanie v obrazovkových pixeloch</translation>
    </message>
    <message>
        <source>Symbology</source>
        <translation>Symbolika</translation>
    </message>
    <message>
        <source>Column 1</source>
        <translation type="obsolete">Stĺpec 1</translation>
    </message>
    <message>
        <source>Table</source>
        <translation>Tabuľka</translation>
    </message>
    <message>
        <source>Add Column</source>
        <translation>Pridať stĺpec</translation>
    </message>
    <message>
        <source>Create / Alter Table</source>
        <translation>Vytvoriť / vymeniť tabuľku</translation>
    </message>
    <message>
        <source>Line width</source>
        <translation>Hrúbka čiary</translation>
    </message>
    <message>
        <source>Marker size</source>
        <translation>Veľkosť značky</translation>
    </message>
    <message>
        <source>Layer</source>
        <translation type="unfinished">Vrstva</translation>
    </message>
    <message>
        <source>Disp</source>
        <translation type="unfinished">Zobraz</translation>
    </message>
    <message>
        <source>Color</source>
        <translation type="unfinished">Farba</translation>
    </message>
    <message>
        <source>Type</source>
        <translation type="unfinished">Typ</translation>
    </message>
    <message>
        <source>Index</source>
        <translation type="unfinished">Index</translation>
    </message>
    <message>
        <source>Column</source>
        <translation type="unfinished">Stĺpec</translation>
    </message>
    <message>
        <source>Length</source>
        <translation type="unfinished">Dĺžka</translation>
    </message>
</context>
<context>
    <name>QgsGrassElementDialog</name>
    <message>
        <source>Cancel</source>
        <translation>Zrušiť</translation>
    </message>
    <message>
        <source>Ok</source>
        <translation>OK</translation>
    </message>
    <message>
        <source>&lt;font color=&apos;red&apos;&gt;Enter a name!&lt;/font&gt;</source>
        <translation type="unfinished">&lt;font color=&apos;red&apos;&gt;Vložte meno!&lt;/font&gt;</translation>
    </message>
    <message>
        <source>&lt;font color=&apos;red&apos;&gt;This is name of the source!&lt;/font&gt;</source>
        <translation type="unfinished">&lt;font color=&apos;red&apos;&gt;Toto je meno zdroja!&lt;/font&gt;</translation>
    </message>
    <message>
        <source>&lt;font color=&apos;red&apos;&gt;Exists!&lt;/font&gt;</source>
        <translation>&lt;font color=&apos;red&apos;&gt;Existuje!&lt;/font&gt;</translation>
    </message>
    <message>
        <source>Overwrite</source>
        <translation type="unfinished">Prepísať</translation>
    </message>
</context>
<context>
    <name>QgsGrassMapcalc</name>
    <message>
        <source>Mapcalc tools</source>
        <translation>Nástroje mapcalc</translation>
    </message>
    <message>
        <source>Add map</source>
        <translation>Pridať mapu</translation>
    </message>
    <message>
        <source>Add constant value</source>
        <translation>Pridať konštantnú hodnotu</translation>
    </message>
    <message>
        <source>Add operator or function</source>
        <translation>Pridať operátor alebo funkciu</translation>
    </message>
    <message>
        <source>Add connection</source>
        <translation>Pridať spojenie</translation>
    </message>
    <message>
        <source>Select item</source>
        <translation>Vybrať položku</translation>
    </message>
    <message>
        <source>Delete selected item</source>
        <translation>Vymazať vybranú položku</translation>
    </message>
    <message>
        <source>Open</source>
        <translation>Otvoriť</translation>
    </message>
    <message>
        <source>Save</source>
        <translation>Uložiť</translation>
    </message>
    <message>
        <source>Save as</source>
        <translation>Uložiť ako</translation>
    </message>
    <message>
        <source>Addition</source>
        <translation>Sčítanie</translation>
    </message>
    <message>
        <source>Subtraction</source>
        <translation>Odčítanie</translation>
    </message>
    <message>
        <source>Multiplication</source>
        <translation>Násobenie</translation>
    </message>
    <message>
        <source>Division</source>
        <translation>Delenie</translation>
    </message>
    <message>
        <source>Modulus</source>
        <translation>Modulo (zvyšok po celočíselnom delení)</translation>
    </message>
    <message>
        <source>Exponentiation</source>
        <translation>Mocnina</translation>
    </message>
    <message>
        <source>Equal</source>
        <translation>Rovnosť</translation>
    </message>
    <message>
        <source>Not equal</source>
        <translation>Nerovnosť</translation>
    </message>
    <message>
        <source>Greater than</source>
        <translation>Väčší než</translation>
    </message>
    <message>
        <source>Greater than or equal</source>
        <translation>Väčší alebo rovný</translation>
    </message>
    <message>
        <source>Less than</source>
        <translation>Menší než</translation>
    </message>
    <message>
        <source>Less than or equal</source>
        <translation>Menší alebo rovný</translation>
    </message>
    <message>
        <source>And</source>
        <translation>A</translation>
    </message>
    <message>
        <source>Or</source>
        <translation>Alebo</translation>
    </message>
    <message>
        <source>Absolute value of x</source>
        <translation>Absolútna hodnota z x</translation>
    </message>
    <message>
        <source>Inverse tangent of x (result is in degrees)</source>
        <translation>Arkus tangens (výsledok v stupňoch)</translation>
    </message>
    <message>
        <source>Inverse tangent of y/x (result is in degrees)</source>
        <translation>Arkus tangens z y/x (výsledok v stupňoch)</translation>
    </message>
    <message>
        <source>Current column of moving window (starts with 1)</source>
        <translation>Aktuálny stĺpec z pohybujúceho sa okna (začína s 1)</translation>
    </message>
    <message>
        <source>Cosine of x (x is in degrees)</source>
        <translation>Kosínus x (x je v stupňoch)</translation>
    </message>
    <message>
        <source>Convert x to double-precision floating point</source>
        <translation>Prekonvertuje x na číslo s pohyblivou rádovou čiarkou s dvojitou presnosťou</translation>
    </message>
    <message>
        <source>Current east-west resolution</source>
        <translation>Aktuálne rozlíšenie vo východo-západnom smere</translation>
    </message>
    <message>
        <source>Exponential function of x</source>
        <translation>Exponenciálna funkcia e na x-tú</translation>
    </message>
    <message>
        <source>x to the power y</source>
        <translation>x na y (mocnina)</translation>
    </message>
    <message>
        <source>Convert x to single-precision floating point</source>
        <translation>Prekonvertuje x na číslo s pohyblivou rádovou čiarkou s jednoduchou presnosťou</translation>
    </message>
    <message>
        <source>Decision: 1 if x not zero, 0 otherwise</source>
        <translation>Rozhodnutie: 1 ak x je nenulové, 0 v ostatných prípadoch</translation>
    </message>
    <message>
        <source>Decision: a if x not zero, 0 otherwise</source>
        <translation>Rozhodnutie: a ak x je nenulové, 0 v ostatných prípadoch</translation>
    </message>
    <message>
        <source>Decision: a if x not zero, b otherwise</source>
        <translation>Rozhodnutie: a ak x je nenulové, b v ostatných prípadoch</translation>
    </message>
    <message>
        <source>Decision: a if x &gt; 0, b if x is zero, c if x &lt; 0</source>
        <translation>Rozhodnutie: a ak x &gt; 0, b ak x je nula, c ak x &lt; 0</translation>
    </message>
    <message>
        <source>Convert x to integer [ truncates ]</source>
        <translation>Prekonvertuje x na celé číslo [ zanedbá desatinnú časť ]</translation>
    </message>
    <message>
        <source>Check if x = NULL</source>
        <translation>Skontroluje či x = NULL (prázdna hodnota)</translation>
    </message>
    <message>
        <source>Natural log of x</source>
        <translation>Prirodzený logaritmus z x</translation>
    </message>
    <message>
        <source>Log of x base b</source>
        <translation>Logaritmus z x pri základe b</translation>
    </message>
    <message>
        <source>Largest value</source>
        <translation>Najväčšia hotnota</translation>
    </message>
    <message>
        <source>Median value</source>
        <translation>Stredná hodnota (medián)</translation>
    </message>
    <message>
        <source>Smallest value</source>
        <translation>Najmenšia hodnota</translation>
    </message>
    <message>
        <source>Mode value</source>
        <translation>Modus</translation>
    </message>
    <message>
        <source>1 if x is zero, 0 otherwise</source>
        <translation>1 ak x je nula, 0 v ostatných pípadoch</translation>
    </message>
    <message>
        <source>Current north-south resolution</source>
        <translation>Aktuálne rozlíšenie v severo-južnom smere</translation>
    </message>
    <message>
        <source>NULL value</source>
        <translation>hodnota NULL (prázdna hodnota)</translation>
    </message>
    <message>
        <source>Random value between a and b</source>
        <translation>Náhodná hodnota medzi a a b</translation>
    </message>
    <message>
        <source>Round x to nearest integer</source>
        <translation>Zaokrúhliť x k najbližšiemu celému číslu</translation>
    </message>
    <message>
        <source>Current row of moving window (Starts with 1)</source>
        <translation>Aktuálny rad pohybujúceho sa okna (začína s 1)</translation>
    </message>
    <message>
        <source>Sine of x (x is in degrees)</source>
        <comment>sin(x)</comment>
        <translation>Sínus x (x je v stupňoch)</translation>
    </message>
    <message>
        <source>Square root of x</source>
        <comment>sqrt(x)</comment>
        <translation>Druhá odmocnina z x</translation>
    </message>
    <message>
        <source>Tangent of x (x is in degrees)</source>
        <comment>tan(x)</comment>
        <translation>Tangens x (x je v stupňoch)</translation>
    </message>
    <message>
        <source>Current x-coordinate of moving window</source>
        <translation>Aktuálna x-ová súradnica pohybujúceho sa okna</translation>
    </message>
    <message>
        <source>Current y-coordinate of moving window</source>
        <translation>Aktuálna y-ová súradnica pohybujúceho sa okna</translation>
    </message>
    <message>
        <source>Warning</source>
        <translation>Upozornenie</translation>
    </message>
    <message>
        <source>Cannot get current region</source>
        <translation>Nemožno zistiť aktuálny región</translation>
    </message>
    <message>
        <source>Cannot check region of map </source>
        <translation type="unfinished">Nemožno skontrolovať región z mapy </translation>
    </message>
    <message>
        <source>Cannot get region of map </source>
        <translation>Nemožno zistiť aktuálny región z mapy </translation>
    </message>
    <message>
        <source>No GRASS raster maps currently in QGIS</source>
        <translation>V QGISe sa momentálne nenachádzajú žiadne rastrové vrstvy GRASSu</translation>
    </message>
    <message>
        <source>Cannot create &apos;mapcalc&apos; directory in current mapset.</source>
        <translation>V aktuálnom súbore máp (mapsete) nemožno vytvoriť adresár &apos;mapcalc&apos;.</translation>
    </message>
    <message>
        <source>New mapcalc</source>
        <translation>Nová schéma pre mapcalc</translation>
    </message>
    <message>
        <source>Enter new mapcalc name:</source>
        <translation type="unfinished">Vložte nové názov schémy pre mapcalc:</translation>
    </message>
    <message>
        <source>Enter vector name</source>
        <translation type="unfinished">Vložte meno vektora</translation>
    </message>
    <message>
        <source>The file already exists. Overwrite? </source>
        <translation>Tento súbor už existuje. Prepísať? </translation>
    </message>
    <message>
        <source>Save mapcalc</source>
        <translation type="unfinished">Uložiť schému mapcalc</translation>
    </message>
    <message>
        <source>File name empty</source>
        <translation>Prázdne meno súboru</translation>
    </message>
    <message>
        <source>Cannot open mapcalc file</source>
        <translation type="unfinished">Nemožno otvoriť súbor mapcalcu</translation>
    </message>
    <message>
        <source>The mapcalc schema (</source>
        <translation type="unfinished">Schéma mapcalcu (</translation>
    </message>
    <message>
        <source>) not found.</source>
        <translation type="unfinished">) nenájdená.</translation>
    </message>
    <message>
        <source>Cannot open mapcalc schema (</source>
        <translation type="unfinished">Nemožno otvoriť schému mapcalcu (</translation>
    </message>
    <message>
        <source>Cannot read mapcalc schema (</source>
        <translation type="unfinished">Nemožno čítať mapcalc schému (</translation>
    </message>
    <message>
        <source>
at line </source>
        <translation type="unfinished">
v riadku </translation>
    </message>
    <message>
        <source> column </source>
        <translation type="unfinished"> stĺpec </translation>
    </message>
    <message>
        <source>Output</source>
        <translation>Výstup</translation>
    </message>
</context>
<context>
    <name>QgsGrassMapcalcBase</name>
    <message>
        <source>Output</source>
        <translation>Výstup</translation>
    </message>
    <message>
        <source>MainWindow</source>
        <translation>MapCalc</translation>
    </message>
</context>
<context>
    <name>QgsGrassModule</name>
    <message>
        <source>Run</source>
        <translation>Spustiť</translation>
    </message>
    <message>
        <source>Stop</source>
        <translation>Zastaviť</translation>
    </message>
    <message>
        <source>Module</source>
        <translation>Modul</translation>
    </message>
    <message>
        <source>Warning</source>
        <translation>Upozornenie</translation>
    </message>
    <message>
        <source>The module file (</source>
        <translation>Súbor modulu (</translation>
    </message>
    <message>
        <source>) not found.</source>
        <translation>) sa nenašiel.</translation>
    </message>
    <message>
        <source>Cannot open module file (</source>
        <translation>Nemožno otvoriť súbor modulu (</translation>
    </message>
    <message>
        <source>)</source>
        <translation>)</translation>
    </message>
    <message>
        <source>Cannot read module file (</source>
        <translation>Nemožno čítať súbor modulu (</translation>
    </message>
    <message>
        <source>):
</source>
        <translation>): 
</translation>
    </message>
    <message>
        <source>
at line </source>
        <translation type="unfinished">
na riadku </translation>
    </message>
    <message>
        <source>Module </source>
        <translation>Modul </translation>
    </message>
    <message>
        <source> not found</source>
        <translation> sa nenašiel</translation>
    </message>
    <message>
        <source>Cannot find man page </source>
        <translation type="unfinished">Nemožno nájsť manuálovú stránku </translation>
    </message>
    <message>
        <source>Not available, cannot open description (</source>
        <translation type="unfinished">Nedostupné, nemožno otvoriť popis (</translation>
    </message>
    <message>
        <source> column </source>
        <translation> stĺpec </translation>
    </message>
    <message>
        <source>Not available, incorrect description (</source>
        <translation type="unfinished">Nedostupné, nesprávny popis (</translation>
    </message>
    <message>
        <source>Cannot get input region</source>
        <translation type="unfinished">Nemožno získať vstupný región</translation>
    </message>
    <message>
        <source>Cannot find module </source>
        <translation>Nemožno nájsť modul </translation>
    </message>
    <message>
        <source>Cannot start module: </source>
        <translation>Nemožno spustiť modul: </translation>
    </message>
    <message>
        <source>&lt;B&gt;Successfully finished&lt;/B&gt;</source>
        <translation>&lt;B&gt;Úspešne dokončené&lt;/B&gt;</translation>
    </message>
    <message>
        <source>&lt;B&gt;Finished with error&lt;/B&gt;</source>
        <translation>&lt;B&gt;Dokončené s chybami&lt;/B&gt;</translation>
    </message>
    <message>
        <source>&lt;B&gt;Module crashed or killed&lt;/B&gt;</source>
        <translation>&lt;B&gt;Modul havaroval alebo bol zabitý&lt;/B&gt;</translation>
    </message>
    <message>
        <source>Use Input Region</source>
        <translation type="unfinished">Použiťregión vstupu</translation>
    </message>
    <message>
        <source>Not available, description not found (</source>
        <translation type="unfinished">Nie je dostupné, popis sa nenašiel (</translation>
    </message>
    <message>
        <source>Please ensure you have the GRASS documentation installed.</source>
        <translation type="unfinished">Prosím uistite sa, že máte nainštalovaná dokumentáciu ku GRASSu.</translation>
    </message>
</context>
<context>
    <name>QgsGrassModuleBase</name>
    <message>
        <source>GRASS Module</source>
        <translation>GRASS - moduly</translation>
    </message>
    <message>
        <source>Options</source>
        <translation>Možnosti</translation>
    </message>
    <message>
        <source>Output</source>
        <translation>Výstup</translation>
    </message>
    <message>
        <source>Run</source>
        <translation>Spustiť</translation>
    </message>
    <message>
        <source>Close</source>
        <translation>Zatvoriť</translation>
    </message>
    <message>
        <source>Manual</source>
        <translation>Manuál</translation>
    </message>
    <message>
        <source>View output</source>
        <translation>Prezrieť výstup</translation>
    </message>
    <message>
        <source>TextLabel</source>
        <translation type="unfinished">TextLabel</translation>
    </message>
</context>
<context>
    <name>QgsGrassModuleField</name>
    <message>
        <source>Attribute field</source>
        <translation type="unfinished">Atribútové pole</translation>
    </message>
</context>
<context>
    <name>QgsGrassModuleFile</name>
    <message>
        <source>File</source>
        <translation>Súbor</translation>
    </message>
    <message>
        <source>:&amp;nbsp;missing value</source>
        <translation>:&amp;nbsp;chýbajúca hodnota</translation>
    </message>
    <message>
        <source>:&amp;nbsp;directory does not exist</source>
        <translation>:&amp;nbsp;adresár neexistuje</translation>
    </message>
</context>
<context>
    <name>QgsGrassModuleGdalInput</name>
    <message>
        <source>Warning</source>
        <translation>Upozornenie</translation>
    </message>
    <message>
        <source>Cannot find layeroption </source>
        <translation type="unfinished">Nemožno nájsť voľbu pre vrstvu </translation>
    </message>
    <message>
        <source>PostGIS driver in OGR does not support schemas!&lt;br&gt;Only the table name will be used.&lt;br&gt;It can result in wrong input if more tables of the same name&lt;br&gt;are present in the database.</source>
        <translation type="unfinished">Ovládač PostGIS nepodporuje schémy!&lt;br&gt;Bude použitý iba názov tabuľky.&lt;br&gt;To môže viesť k chybnému výsledku pokiaľ má viacero tabuliek to isté meno&lt;br&gt;v jednej databáze.</translation>
    </message>
    <message>
        <source>:&amp;nbsp;no input</source>
        <translation>:&amp;nbsp;žiadny vstup</translation>
    </message>
    <message>
        <source>Cannot find whereoption </source>
        <translation type="unfinished">Nemožno nájsť nastavenie where </translation>
    </message>
</context>
<context>
    <name>QgsGrassModuleInput</name>
    <message>
        <source>Warning</source>
        <translation>Upozornenie</translation>
    </message>
    <message>
        <source>Cannot find typeoption </source>
        <translation type="unfinished">Nemožno nájsť voľbu type </translation>
    </message>
    <message>
        <source>Cannot find values for typeoption </source>
        <translation type="unfinished">Nemožno nájsť hodnoty pre voľbu typ </translation>
    </message>
    <message>
        <source>Cannot find layeroption </source>
        <translation type="unfinished">Nemožno nájsť voľbu vrstva </translation>
    </message>
    <message>
        <source>GRASS element </source>
        <translation type="unfinished">GRASS prvok </translation>
    </message>
    <message>
        <source> not supported</source>
        <translation type="unfinished"> nie je podporovaný</translation>
    </message>
    <message>
        <source>Use region of this map</source>
        <translation type="unfinished">Použite región tejto mapy</translation>
    </message>
    <message>
        <source>:&amp;nbsp;no input</source>
        <translation>:&amp;nbsp;žiaden vstup</translation>
    </message>
</context>
<context>
    <name>QgsGrassModuleOption</name>
    <message>
        <source>:&amp;nbsp;missing value</source>
        <translation>:&amp;nbsp;chýbajúca hodnota</translation>
    </message>
</context>
<context>
    <name>QgsGrassModuleSelection</name>
    <message>
        <source>Attribute field</source>
        <translation type="unfinished">Atribútové pole</translation>
    </message>
</context>
<context>
    <name>QgsGrassModuleStandardOptions</name>
    <message>
        <source>Warning</source>
        <translation>Upozornenie</translation>
    </message>
    <message>
        <source>Cannot find module </source>
        <translation type="unfinished">Nemožno nájsť modul </translation>
    </message>
    <message>
        <source>Cannot start module </source>
        <translation>Nemožno spustiť modul </translation>
    </message>
    <message>
        <source>Cannot read module description (</source>
        <translation>Nemožno prečítať popis k modulu (</translation>
    </message>
    <message>
        <source>):
</source>
        <translation>):
</translation>
    </message>
    <message>
        <source>
at line </source>
        <translation type="unfinished">
na riadku </translation>
    </message>
    <message>
        <source> column </source>
        <translation type="unfinished"> stĺpec </translation>
    </message>
    <message>
        <source>Cannot find key </source>
        <translation type="unfinished">Nemožno nájsť kľúč </translation>
    </message>
    <message>
        <source>Item with id </source>
        <translation type="unfinished">Položka s id </translation>
    </message>
    <message>
        <source> not found</source>
        <translation type="unfinished"> sa nenašla</translation>
    </message>
    <message>
        <source>Cannot get current region</source>
        <translation type="unfinished">Nemožno zistiť aktuálny región</translation>
    </message>
    <message>
        <source>Cannot check region of map </source>
        <translation type="unfinished">Nemožno skontrolovať región mapy </translation>
    </message>
    <message>
        <source>Cannot set region of map </source>
        <translation type="unfinished">Nemožno nastaviť región mapy </translation>
    </message>
</context>
<context>
    <name>QgsGrassNewMapset</name>
    <message>
        <source>GRASS database</source>
        <translation type="obsolete">Databáza GRASSu</translation>
    </message>
    <message>
        <source>GRASS location</source>
        <translation type="obsolete">Lokalita GRASSu</translation>
    </message>
    <message>
        <source>Projection</source>
        <translation type="obsolete">Mapové zobrazenie</translation>
    </message>
    <message>
        <source>Default GRASS Region</source>
        <translation type="obsolete">Predvolený región GRASSu</translation>
    </message>
    <message>
        <source>Mapset</source>
        <translation type="obsolete">Súbor máp (mapset)</translation>
    </message>
    <message>
        <source>Create New Mapset</source>
        <translation type="obsolete">Vytvoriť nový súbor máp (mapset)</translation>
    </message>
    <message>
        <source>Tree</source>
        <translation type="obsolete">Strom</translation>
    </message>
    <message>
        <source>Comment</source>
        <translation type="obsolete">Komentár</translation>
    </message>
    <message>
        <source>Database</source>
        <translation>Databáza</translation>
    </message>
    <message>
        <source>Location 2</source>
        <translation>Lokalita 2</translation>
    </message>
    <message>
        <source>User&apos;s mapset</source>
        <translation type="unfinished">Uživateľský mapset</translation>
    </message>
    <message>
        <source>System mapset</source>
        <translation type="unfinished">Systémový mapset</translation>
    </message>
    <message>
        <source>Location 1</source>
        <translation>Lokalita 1</translation>
    </message>
    <message>
        <source>Owner</source>
        <translation type="obsolete">Vlastník</translation>
    </message>
    <message>
        <source>Enter path to GRASS database</source>
        <translation>Zadajte cestu k databáze GRASSu</translation>
    </message>
    <message>
        <source>The directory doesn&apos;t exist!</source>
        <translation type="unfinished">Tento adresár neexistuje!</translation>
    </message>
    <message>
        <source>No writable locations, the database not writable!</source>
        <translation type="unfinished">Žiadne zapisovateľné lokality, do databázy nemožno zapisovať!</translation>
    </message>
    <message>
        <source>Enter location name!</source>
        <translation type="unfinished">Zadajte názov lokality!</translation>
    </message>
    <message>
        <source>The location exists!</source>
        <translation type="unfinished">Táto lokalita existuje!</translation>
    </message>
    <message>
        <source>Selected projection is not supported by GRASS!</source>
        <translation type="unfinished">Vybrané mapové zobraznie nie je podporované GRASSom!</translation>
    </message>
    <message>
        <source>Warning</source>
        <translation>Upozornenie</translation>
    </message>
    <message>
        <source>Cannot create projection.</source>
        <translation type="unfinished">Nemožno vytvoriť mapové zobrazenie.</translation>
    </message>
    <message>
        <source>Cannot reproject previously set region, default region set.</source>
        <translation type="unfinished">Nemožno prezobraziť predtým nastavený región, predvolený región bol nastavený</translation>
    </message>
    <message>
        <source>North must be greater than south</source>
        <translation type="unfinished">Sever musí mať väčšiu hodnotu ako juh</translation>
    </message>
    <message>
        <source>East must be greater than west</source>
        <translation type="unfinished">Východ musí mať väčšiu hodnotu ako západ</translation>
    </message>
    <message>
        <source>Regions file (</source>
        <translation type="unfinished">Súbor s regiónom (</translation>
    </message>
    <message>
        <source>) not found.</source>
        <translation type="unfinished">) sa nenašiel.</translation>
    </message>
    <message>
        <source>Cannot open locations file (</source>
        <translation type="unfinished">Nemožno otvoriť súbor s lokalitami (</translation>
    </message>
    <message>
        <source>)</source>
        <translation>)</translation>
    </message>
    <message>
        <source>Cannot read locations file (</source>
        <translation type="unfinished">Nemožno čítať súbor s lokalitami (</translation>
    </message>
    <message>
        <source>):
</source>
        <translation type="unfinished">):
</translation>
    </message>
    <message>
        <source>
at line </source>
        <translation type="unfinished">
v riadku </translation>
    </message>
    <message>
        <source> column </source>
        <translation type="unfinished"> stĺpec </translation>
    </message>
    <message>
        <source>Cannot create QgsSpatialRefSys</source>
        <translation type="obsolete">Nemožno vytvoriť QgsSpatialRefSys</translation>
    </message>
    <message>
        <source>Cannot reproject selected region.</source>
        <translation type="unfinished">Nemožno reprojektovať vybraný región.</translation>
    </message>
    <message>
        <source>Cannot reproject region</source>
        <translation type="unfinished">Nemožno prezobraziť región</translation>
    </message>
    <message>
        <source>Enter mapset name.</source>
        <translation type="unfinished">Zadajte názov súboru máp (mapsetu)</translation>
    </message>
    <message>
        <source>The mapset already exists</source>
        <translation type="unfinished">Súbor máp (mapset) už existuje</translation>
    </message>
    <message>
        <source>Database: </source>
        <translation>Databáza: </translation>
    </message>
    <message>
        <source>Location: </source>
        <translation>Lokalita: </translation>
    </message>
    <message>
        <source>Mapset: </source>
        <translation type="unfinished">Súbor máp (mapset): </translation>
    </message>
    <message>
        <source>Create location</source>
        <translation type="unfinished">Vytvoriť lokalitu</translation>
    </message>
    <message>
        <source>Cannot create new location: </source>
        <translation type="unfinished">Nemožno vytvoriť novú lokalitu: </translation>
    </message>
    <message>
        <source>Create mapset</source>
        <translation type="unfinished">Vytvoriť mapset</translation>
    </message>
    <message>
        <source>Cannot open DEFAULT_WIND</source>
        <translation type="unfinished">Nemožno otvoriť DEFAULT_WIND</translation>
    </message>
    <message>
        <source>Cannot open WIND</source>
        <translation type="unfinished">Nemožno otvoriť WIND</translation>
    </message>
    <message>
        <source>New mapset</source>
        <translation type="unfinished">Nový súbor máp (mapset)</translation>
    </message>
    <message>
        <source>New mapset successfully created, but cannot be opened: </source>
        <translation type="unfinished">Nový mapset bol úspečšne vytvorený ale nemožno ho otvoriť: </translation>
    </message>
    <message>
        <source>New mapset successfully created and set as current working mapset.</source>
        <translation type="unfinished">Nový mapset bol úspešne vytvorený a nastavený ako súčasný pracovný mapset.</translation>
    </message>
    <message>
        <source>Cannot create new mapset directory</source>
        <translation type="unfinished">Nemožno vytvoriť nový adresár pre mapset (zbierku máp)</translation>
    </message>
    <message>
        <source>Cannot create QgsCoordinateReferenceSystem</source>
        <translation type="unfinished">Nemožno vytvoriť QgsCoordinateReferenceSystem</translation>
    </message>
</context>
<context>
    <name>QgsGrassNewMapsetBase</name>
    <message>
        <source>Select existing directory or create a new one:</source>
        <translation>Vyberte existujúci adresár alebo vytvorte nový:</translation>
    </message>
    <message>
        <source>Database:</source>
        <translation>Databáza:</translation>
    </message>
    <message>
        <source>...</source>
        <translation type="obsolete">...</translation>
    </message>
    <message>
        <source>Example directory tree:</source>
        <translation>Príklad adresárovej štruktúry:</translation>
    </message>
    <message>
        <source>Column 1</source>
        <translation type="obsolete">Stĺpec 1</translation>
    </message>
    <message>
        <source>Database Error</source>
        <translation>Chyba databázy</translation>
    </message>
    <message>
        <source>Location</source>
        <translation>Lokalita</translation>
    </message>
    <message>
        <source>Select location</source>
        <translation>Vybrať lokalitu</translation>
    </message>
    <message>
        <source>Create new location</source>
        <translation>Vytvoriť novú lokalitu</translation>
    </message>
    <message>
        <source>Location Error</source>
        <translation>Chyba v lokalite</translation>
    </message>
    <message>
        <source>Coordinate system</source>
        <translation>Súradnicový systém</translation>
    </message>
    <message>
        <source>Not defined</source>
        <translation>Nedefinované</translation>
    </message>
    <message>
        <source>Projection</source>
        <translation>Mapové zobrazenie</translation>
    </message>
    <message>
        <source>Projection Error</source>
        <translation>Chyba mapového zobrazenia</translation>
    </message>
    <message>
        <source>N</source>
        <translation>S</translation>
    </message>
    <message>
        <source>W</source>
        <translation>Z</translation>
    </message>
    <message>
        <source>E</source>
        <translation>V</translation>
    </message>
    <message>
        <source>S</source>
        <translation>J</translation>
    </message>
    <message>
        <source>Region Error</source>
        <translation>Chyba v regióne</translation>
    </message>
    <message>
        <source>Set current QGIS extent</source>
        <translation>Nastaviť na aktuálny rozsah v QGISe</translation>
    </message>
    <message>
        <source>Set</source>
        <translation>Nastaviť</translation>
    </message>
    <message>
        <source>New mapset:</source>
        <translation>Nový súbor máp (mapset):</translation>
    </message>
    <message>
        <source>Mapset Error</source>
        <translation>Chyba v súbore máp</translation>
    </message>
    <message>
        <source>&lt;p align=&quot;center&quot;&gt;Existing masets&lt;/p&gt;</source>
        <translation>&lt;p align=&quot;center&quot;&gt;Existujúce súbory máp (mapsety)&lt;/p&gt;</translation>
    </message>
    <message>
        <source>Location:</source>
        <translation>Lokalita:</translation>
    </message>
    <message>
        <source>Mapset:</source>
        <translation>Súbor máp:</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;/head&gt;&lt;body style=&quot; white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;The GRASS location is a collection of maps for a particular territory or project.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;/head&gt;&lt;body style=&quot; white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Lokalita v GRASSe je chápaná zbierka máp pre určité územie alebo projekt.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;/head&gt;&lt;body style=&quot; white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;The GRASS mapset is a collection of maps used by one user. &lt;/p&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;A user can read maps from all mapsets in the location but &lt;/p&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;he can open for writing only his mapset (owned by user).&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;/head&gt;&lt;body style=&quot; white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Súbor máp GRASS-u (GRASS mapset) je zbierka máp používaných jedným užívateľom. &lt;/p&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Užívateľ môže čítať mapy zo všetkých mapových súprav v lokalite, ale &lt;/p&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;otvoriť zapisovovať môže len do svojho súboru máp (musí byť vlastníkom daného súboru).&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;/head&gt;&lt;body style=&quot; white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;GRASS data are stored in tree directory structure.&lt;/p&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;The GRASS database is the top-level directory in this tree structure.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;/head&gt;&lt;body style=&quot; white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Údaje GRASSu sú uložené v stromovej adresárovej štruktúre.&lt;/p&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Databáza GRASSu je adresár najvyššej úrovne tejto stromovej adresárovej štruktúry.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;/head&gt;&lt;body style=&quot; white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;The GRASS region defines a workspace for raster modules. The default region is valid for one location. It is possible to set a different region in each mapset. &lt;/p&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;It is possible to change the default location region later.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;/head&gt;&lt;body style=&quot; white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Región GRASSu určuje pracovný priestor pre rastrové moduly. Predvolený región je platný pre jednu lokalitu. Pre každý súbor máp (mapset) je možné nastaviť iný región &lt;/p&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Predvolený región pre lokalitu je možné neskôr zmeniť.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;GRASS data are stored in tree directory structure.&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;The GRASS database is the top-level directory in this tree structure.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;Údaje GRASSu sú uložené v stromovej adresárovej štruktúre.&lt;/p&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;Databáza GRASSu je adresár najvyššej úrovne tejto stromovej adresárovej štruktúry.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;The GRASS location is a collection of maps for a particular territory or project.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;Lokalita v GRASSe je chápaná zbierka máp pre určité územie alebo projekt.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;The GRASS region defines a workspace for raster modules. The default region is valid for one location. It is possible to set a different region in each mapset. &lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;It is possible to change the default location region later.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;Región GRASSu určuje pracovný priestor pre rastrové moduly. Predvolený región je platný pre jednu lokalitu. Pre každý súbor máp (mapset) je možné nastaviť iný región.&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;Predvolený región pre lokalitu je možné neskôr zmeniť.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;The GRASS mapset is a collection of maps used by one user. &lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;A user can read maps from all mapsets in the location but &lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;he can open for writing only his mapset (owned by user).&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;Súbor máp GRASS-u (GRASS mapset) je zbierka máp používaných jedným užívateľom. &lt;/p&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;Užívateľ môže čítať mapy zo všetkých mapových súprav v lokalite, ale &lt;/p&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;zapisovovať môže len do svojho súboru máp (musí byť vlastníkom daného súboru).&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>New Mapset</source>
        <translation type="unfinished">Nový súbor máp (mapset)</translation>
    </message>
    <message>
        <source>GRASS Database</source>
        <translation type="unfinished">Databáza GRASSu</translation>
    </message>
    <message>
        <source>Tree</source>
        <translation type="unfinished">Strom</translation>
    </message>
    <message>
        <source>Comment</source>
        <translation type="unfinished">Komentár</translation>
    </message>
    <message>
        <source>&lt;!DOCTYPE HTML PUBLIC &quot;-//W3C//DTD HTML 4.0//EN&quot; &quot;http://www.w3.org/TR/REC-html40/strict.dtd&quot;&gt;
&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Lucida Grande&apos;; font-size:13pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;GRASS data are stored in tree directory structure. The GRASS database is the top-level directory in this tree structure.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="unfinished">&lt;!DOCTYPE HTML PUBLIC &quot;-//W3C//DTD HTML 4.0//EN&quot; &quot;http://www.w3.org/TR/REC-html40/strict.dtd&quot;&gt;
&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Lucida Grande&apos;; font-size:13pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;Údaje GRASSu sú uložené v stromovej adresárovej štruktúre. Databáza GRASSu je adresár najvyššej úrovne tejto stromovej adresárovej štruktúry.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>Browse...</source>
        <translation>Prechádzať...</translation>
    </message>
    <message>
        <source>GRASS Location</source>
        <translation type="unfinished">Lokalita GRASSu</translation>
    </message>
    <message>
        <source>&lt;!DOCTYPE HTML PUBLIC &quot;-//W3C//DTD HTML 4.0//EN&quot; &quot;http://www.w3.org/TR/REC-html40/strict.dtd&quot;&gt;
&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Lucida Grande&apos;; font-size:13pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;The GRASS location is a collection of maps for a particular territory or project.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="unfinished">&lt;!DOCTYPE HTML PUBLIC &quot;-//W3C//DTD HTML 4.0//EN&quot; &quot;http://www.w3.org/TR/REC-html40/strict.dtd&quot;&gt;
&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Lucida Grande&apos;; font-size:13pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;Lokalita v GRASSe je chápaná zbierka máp pre určité územie alebo projekt.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>Default GRASS Region</source>
        <translation>Predvolený región GRASSu</translation>
    </message>
    <message>
        <source>&lt;!DOCTYPE HTML PUBLIC &quot;-//W3C//DTD HTML 4.0//EN&quot; &quot;http://www.w3.org/TR/REC-html40/strict.dtd&quot;&gt;
&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Lucida Grande&apos;; font-size:13pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;The GRASS region defines a workspace for raster modules. The default region is valid for one location. It is possible to set a different region in each mapset. It is possible to change the default location region later.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="unfinished">&lt;!DOCTYPE HTML PUBLIC &quot;-//W3C//DTD HTML 4.0//EN&quot; &quot;http://www.w3.org/TR/REC-html40/strict.dtd&quot;&gt;
&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Lucida Grande&apos;; font-size:13pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;Región GRASSu určuje pracovný priestor pre rastrové moduly. Predvolený región je platný pre jednu lokalitu. Pre každý súbor máp (mapset) je možné nastaviť iný región. It is possible to change the default location region later. Predvolený región pre lokalitu je možné neskôr zmeniť.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>Mapset</source>
        <translation type="unfinished">Súbor máp (mapset)</translation>
    </message>
    <message>
        <source>Owner</source>
        <translation>Vlastník</translation>
    </message>
    <message>
        <source>&lt;!DOCTYPE HTML PUBLIC &quot;-//W3C//DTD HTML 4.0//EN&quot; &quot;http://www.w3.org/TR/REC-html40/strict.dtd&quot;&gt;
&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Lucida Grande&apos;; font-size:13pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;The GRASS mapset is a collection of maps used by one user. A user can read maps from all mapsets in the location but he can open for writing only his mapset (owned by user).&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="unfinished">&lt;!DOCTYPE HTML PUBLIC &quot;-//W3C//DTD HTML 4.0//EN&quot; &quot;http://www.w3.org/TR/REC-html40/strict.dtd&quot;&gt;&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;p, li { white-space: pre-wrap; }&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Lucida Grande&apos;; font-size:13pt; font-weight:400; font-style:normal;&quot;&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;Súbor máp GRASS-u (GRASS mapset) je zbierka máp používaných jedným užívateľom. Užívateľ môže čítať mapy zo všetkých mapových súprav v lokalite, ale zapisovovať môže len do svojho súboru máp (musí byť vlastníkom daného súboru).&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>Create New Mapset</source>
        <translation type="unfinished">Vytvoriť nový súbor máp (mapset)</translation>
    </message>
</context>
<context>
    <name>QgsGrassPlugin</name>
    <message>
        <source>GRASS</source>
        <translation>GRASS</translation>
    </message>
    <message>
        <source>&amp;GRASS</source>
        <translation>&amp;GRASS</translation>
    </message>
    <message>
        <source>Open mapset</source>
        <translation>Otvoriť súbor máp (mapset)</translation>
    </message>
    <message>
        <source>New mapset</source>
        <translation>Nový súbor máp (mapset)</translation>
    </message>
    <message>
        <source>Close mapset</source>
        <translation>Zatvoriť súbor máp (mapset)</translation>
    </message>
    <message>
        <source>Add GRASS vector layer</source>
        <translation>Pridať vektorovú vrstvu GRASSu</translation>
    </message>
    <message>
        <source>Add GRASS raster layer</source>
        <translation>Pridať rastrovú vrstvu GRASSu</translation>
    </message>
    <message>
        <source>Open GRASS tools</source>
        <translation>Otvoriť nástroje GRASSu</translation>
    </message>
    <message>
        <source>Display Current Grass Region</source>
        <translation>Zobraziť aktuálny región GRASSu</translation>
    </message>
    <message>
        <source>Edit Current Grass Region</source>
        <translation>Upraviť aktuálny región GRASSu</translation>
    </message>
    <message>
        <source>Edit Grass Vector layer</source>
        <translation>Upraviť vektorovú vrstvu GRASSu</translation>
    </message>
    <message>
        <source>Adds a GRASS vector layer to the map canvas</source>
        <translation>Na mapové plátno pridá vektorovú vrstvu GRASSu</translation>
    </message>
    <message>
        <source>Adds a GRASS raster layer to the map canvas</source>
        <translation>Na mapové plátno pridá rastrovú vrstvu GRASSu</translation>
    </message>
    <message>
        <source>Displays the current GRASS region as a rectangle on the map canvas</source>
        <translation>Zobrazí na mapovom plátne aktuálny región GRASSu ako obdĺžnik</translation>
    </message>
    <message>
        <source>Edit the current GRASS region</source>
        <translation>Upraví aktuálny región GRASSu</translation>
    </message>
    <message>
        <source>Edit the currently selected GRASS vector layer.</source>
        <translation>Upraví vybranú vektorovú vrstvu GRASSu.</translation>
    </message>
    <message>
        <source>GrassVector</source>
        <translation type="unfinished">GrassVektor</translation>
    </message>
    <message>
        <source>0.1</source>
        <translation>0.1</translation>
    </message>
    <message>
        <source>GRASS layer</source>
        <translation>Vrstva GRASSu</translation>
    </message>
    <message>
        <source>Warning</source>
        <translation>Upozornenie</translation>
    </message>
    <message>
        <source>GRASS Edit is already running.</source>
        <translation type="unfinished">GRASS úpravy už prebiehajú.</translation>
    </message>
    <message>
        <source>Cannot create new vector: </source>
        <translation type="unfinished">Nemožno vytvoriť nový vektor: </translation>
    </message>
    <message>
        <source>New vector created but cannot be opened by data provider.</source>
        <translation type="unfinished">Nový vektor vytvorený, ale nemôže byť otvorený správcom údajov.</translation>
    </message>
    <message>
        <source>Cannot start editing.</source>
        <translation type="unfinished">Nemožno začať upravovať.</translation>
    </message>
    <message>
        <source>GISDBASE, LOCATION_NAME or MAPSET is not set, cannot display current region.</source>
        <translation type="unfinished">GISDBASE, LOCATION_NAME alebo MAPSET nie je nastavený, nemožno zobraziť aktuálny región.</translation>
    </message>
    <message>
        <source>Cannot read current region: </source>
        <translation type="unfinished">Nemožno čítať aktuálny región: </translation>
    </message>
    <message>
        <source>Cannot open the mapset. </source>
        <translation type="unfinished">Nemožno otvoriť tento súbopr máp. </translation>
    </message>
    <message>
        <source>Cannot close mapset. </source>
        <translation type="unfinished">Nemožno zavrieť súbor máp. </translation>
    </message>
    <message>
        <source>Cannot close current mapset. </source>
        <translation type="unfinished">Nemožno zavrieť aktuálny súbor máp. </translation>
    </message>
    <message>
        <source>Cannot open GRASS mapset. </source>
        <translation type="unfinished">Nemožno otvoriť súbor máp GRASSu. </translation>
    </message>
    <message>
        <source>Create new Grass Vector</source>
        <translation type="unfinished">Vytvoriť novú vektor. vrstvu GRASSu</translation>
    </message>
    <message>
        <source>New vector name</source>
        <translation type="unfinished">Nové meno vektora</translation>
    </message>
</context>
<context>
    <name>QgsGrassRegion</name>
    <message>
        <source>Warning</source>
        <translation>Upozornenie</translation>
    </message>
    <message>
        <source>GISDBASE, LOCATION_NAME or MAPSET is not set, cannot display current region.</source>
        <translation type="unfinished">GISDBASE, LOCATION_NAME alebo MAPSET nie je nastavený, nemožno zobraziť aktuálny región.</translation>
    </message>
    <message>
        <source>Cannot read current region: </source>
        <translation type="unfinished">Nemožno čítať aktuálny región: </translation>
    </message>
    <message>
        <source>Cannot write region</source>
        <translation type="unfinished">Nemožno zapísať región</translation>
    </message>
</context>
<context>
    <name>QgsGrassRegionBase</name>
    <message>
        <source>GRASS Region Settings</source>
        <translation>Nastavenie pracovnej oblasti GRASSu</translation>
    </message>
    <message>
        <source>N</source>
        <translation>S</translation>
    </message>
    <message>
        <source>W</source>
        <translation>Z</translation>
    </message>
    <message>
        <source>E</source>
        <translation>V</translation>
    </message>
    <message>
        <source>S</source>
        <translation>J</translation>
    </message>
    <message>
        <source>N-S Res</source>
        <translation>Rozlíšenie S-J</translation>
    </message>
    <message>
        <source>Rows</source>
        <translation>Riadky</translation>
    </message>
    <message>
        <source>Cols</source>
        <translation>Stĺpce</translation>
    </message>
    <message>
        <source>E-W Res</source>
        <translation>Rozlíšenie V-Z</translation>
    </message>
    <message>
        <source>Color</source>
        <translation>Farba</translation>
    </message>
    <message>
        <source>Width</source>
        <translation>Hrúbka</translation>
    </message>
    <message>
        <source>OK</source>
        <translation>OK</translation>
    </message>
    <message>
        <source>Cancel</source>
        <translation>Zrušiť</translation>
    </message>
</context>
<context>
    <name>QgsGrassSelect</name>
    <message>
        <source>Select GRASS Vector Layer</source>
        <translation>Vyberte vektorovú vrstvu GRASSu</translation>
    </message>
    <message>
        <source>Select GRASS Raster Layer</source>
        <translation>Vyberte rastrovú vrstvu GRASSu</translation>
    </message>
    <message>
        <source>Select GRASS mapcalc schema</source>
        <translation>Vyberte schému mapcalcu</translation>
    </message>
    <message>
        <source>Select GRASS Mapset</source>
        <translation>Vyberte súbor máp (mapset) GRASSu</translation>
    </message>
    <message>
        <source>Warning</source>
        <translation>Upozornenie</translation>
    </message>
    <message>
        <source>Cannot open vector on level 2 (topology not available).</source>
        <translation>Nemožno otvoriť vektor na úrovni 2 (nie je dostupná topológia).</translation>
    </message>
    <message>
        <source>Choose existing GISDBASE</source>
        <translation>Vyberte existujúcu GISDBASE</translation>
    </message>
    <message>
        <source>Wrong GISDBASE, no locations available.</source>
        <translation>Chybná GISDBASE, nie sú dostupné žiadne lokality (locations).</translation>
    </message>
    <message>
        <source>Wrong GISDBASE</source>
        <translation>Chybná GISDBASE</translation>
    </message>
    <message>
        <source>Select a map.</source>
        <translation>Vyberte mapu.</translation>
    </message>
    <message>
        <source>No map</source>
        <translation>Žiadna mapa</translation>
    </message>
    <message>
        <source>No layer</source>
        <translation>Žiadna vrstva</translation>
    </message>
    <message>
        <source>No layers available in this map</source>
        <translation>V tejto mape sa nenachádzajú žiadne vrstvy</translation>
    </message>
</context>
<context>
    <name>QgsGrassSelectBase</name>
    <message>
        <source>Gisdbase</source>
        <translation>Gisdbase</translation>
    </message>
    <message>
        <source>Location</source>
        <translation>Lokalita</translation>
    </message>
    <message>
        <source>Browse</source>
        <translation>Prechádzať</translation>
    </message>
    <message>
        <source>Mapset</source>
        <translation>Mapset</translation>
    </message>
    <message>
        <source>Map name</source>
        <translation>Meno mapy</translation>
    </message>
    <message>
        <source>Layer</source>
        <translation>Vrstva</translation>
    </message>
    <message>
        <source>OK</source>
        <translation>OK</translation>
    </message>
    <message>
        <source>Select or type map name (wildcards &apos;*&apos; and &apos;?&apos; accepted for rasters)</source>
        <translation>Vyberte alebo napíšte meno mapy (pre rastre budú akceptované aj označovacie konvencie &apos;*&apos; a &apos;?&apos;)</translation>
    </message>
    <message>
        <source>Add GRASS Layer</source>
        <translation>Pridať vrstvu GRASSu</translation>
    </message>
    <message>
        <source>Cancel</source>
        <translation>Zrušiť</translation>
    </message>
</context>
<context>
    <name>QgsGrassShellBase</name>
    <message>
        <source>GRASS Shell</source>
        <translation>Príkazový riadok GRASS-u</translation>
    </message>
    <message>
        <source>Close</source>
        <translation>Zatvoriť</translation>
    </message>
</context>
<context>
    <name>QgsGrassTools</name>
    <message>
        <source>Modules</source>
        <translation type="obsolete">Moduly</translation>
    </message>
    <message>
        <source>Browser</source>
        <translation>Prehliadač</translation>
    </message>
    <message>
        <source>GRASS Tools</source>
        <translation>Nástroje GRASSu</translation>
    </message>
    <message>
        <source>GRASS Tools: </source>
        <translation>Nástroje GRASSu: </translation>
    </message>
    <message>
        <source>Warning</source>
        <translation>Upozornenie</translation>
    </message>
    <message>
        <source>Cannot find MSYS (</source>
        <translation type="obsolete">Nemožno nájsť MSYS (</translation>
    </message>
    <message>
        <source>GRASS Shell is not compiled.</source>
        <translation type="unfinished">GRASS Shell nie je skompilovaný.</translation>
    </message>
    <message>
        <source>The config file (</source>
        <translation>Konfiguračný súbor (</translation>
    </message>
    <message>
        <source>) not found.</source>
        <translation>) sa nenašiel.</translation>
    </message>
    <message>
        <source>Cannot open config file (</source>
        <translation type="unfinished">Nemožno otvoriť konfiguračný súbor (</translation>
    </message>
    <message>
        <source>)</source>
        <translation>)</translation>
    </message>
    <message>
        <source>Cannot read config file (</source>
        <translation>Nemožno čítať konfiguračný súbor (</translation>
    </message>
    <message>
        <source>
at line </source>
        <translation type="unfinished">
na riadku </translation>
    </message>
    <message>
        <source> column </source>
        <translation type="unfinished"> stĺpec </translation>
    </message>
    <message>
        <source>Cannot start command shell (%1)</source>
        <translation type="unfinished">Nemožno spustiť príkazový riadok (%1)</translation>
    </message>
</context>
<context>
    <name>QgsGrassToolsBase</name>
    <message>
        <source>Grass Tools</source>
        <translation>Nástroje GRASSu</translation>
    </message>
    <message>
        <source>Modules Tree</source>
        <translation>Strom modulov</translation>
    </message>
    <message>
        <source>1</source>
        <translation type="unfinished">1</translation>
    </message>
    <message>
        <source>Modules List</source>
        <translation>Zoznam modulov</translation>
    </message>
</context>
<context>
    <name>QgsGridMakerPlugin</name>
    <message>
        <source>&amp;Graticules</source>
        <translation>&amp;Súradnicové siete</translation>
    </message>
    <message>
        <source>Creates a graticule (grid) and stores the result as a shapefile</source>
        <translation>Vytvorí súradnicovú sieť (mriežku) a výsledok uloží ako súbor shape</translation>
    </message>
    <message>
        <source>&amp;Graticule Creator</source>
        <translation>&amp;Tvorba súradnicovej siete</translation>
    </message>
</context>
<context>
    <name>QgsGridMakerPluginGui</name>
    <message>
        <source>QGIS - Grid Maker</source>
        <translation type="unfinished">QGIS - Tvorba siete</translation>
    </message>
    <message>
        <source>Choose a filename to save under</source>
        <translation type="obsolete">Vyberte meno súboru do ktorého sa bude ukladať</translation>
    </message>
    <message>
        <source>ESRI Shapefile (*.shp)</source>
        <translation>Súbory ESRI Shape (*.shp)</translation>
    </message>
    <message>
        <source>Please enter the file name before pressing OK!</source>
        <translation type="unfinished">Pred kliknutím na OK, zadajte meno súboru!</translation>
    </message>
    <message>
        <source>Please enter intervals before pressing OK!</source>
        <translation type="unfinished">Prosím zadajte intervaly pred stlačením OK!</translation>
    </message>
    <message>
        <source>Choose a file name to save under</source>
        <translation type="unfinished">Vyberte meno súboru do ktorého sa bude ukladať</translation>
    </message>
</context>
<context>
    <name>QgsGridMakerPluginGuiBase</name>
    <message>
        <source>Graticule Builder</source>
        <translation>Tvorba siete zemepisných súradníc</translation>
    </message>
    <message>
        <source>Type</source>
        <translation>Typ</translation>
    </message>
    <message>
        <source>Point</source>
        <translation>Bod</translation>
    </message>
    <message>
        <source>Polygon</source>
        <translation>Polygón</translation>
    </message>
    <message>
        <source>Origin (lower left)</source>
        <translation>Počiatok (ľavý dolný roh)</translation>
    </message>
    <message>
        <source>End point (upper right)</source>
        <translation>Koncový bod (pravý horný roh)</translation>
    </message>
    <message>
        <source>Output (shape) file</source>
        <translation>Výstupný (shape) súbor</translation>
    </message>
    <message>
        <source>Save As...</source>
        <translation>Uložiť ako...</translation>
    </message>
    <message>
        <source>QGIS Graticule Creator</source>
        <translation type="unfinished">QGIS Tvorba súradnicovej siete</translation>
    </message>
    <message>
        <source>Graticle size</source>
        <translation type="unfinished">Veľkosť siete</translation>
    </message>
    <message>
        <source>Y Interval:</source>
        <translation type="unfinished">Interval v smere Y:</translation>
    </message>
    <message>
        <source>X Interval:</source>
        <translation type="unfinished">Interval v smere X:</translation>
    </message>
    <message>
        <source>Y</source>
        <translation>Y</translation>
    </message>
    <message>
        <source>X</source>
        <translation>X</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;
&lt;p style=&quot; margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Arial&apos;; font-size:11pt;&quot;&gt;&lt;span style=&quot; font-size:10pt;&quot;&gt;This plugin will help you to build a graticule shapefile that you can use as an overlay within your qgis map viewer.&lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Arial&apos;; font-size:10pt;&quot;&gt;Please enter all units in decimal degrees&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;
&lt;p style=&quot; margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Arial&apos;; font-size:11pt;&quot;&gt;&lt;span style=&quot; font-size:10pt;&quot;&gt;Tento zásuvný modul vám pomôže zostrojiť súbor shape so zemepisnou sieťou, ktorú možno preložiť cez mapy vo vašom mapovom prehliadači QGIS.&lt;/span&gt;&lt;/p&gt;&lt;p style=&quot; margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Arial&apos;; font-size:10pt;&quot;&gt;Všetky jednotky zadávajte v desiatkových stupňoch&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Arial&apos;; font-size:11pt;&quot;&gt;&lt;span style=&quot; font-size:10pt;&quot;&gt;This plugin will help you to build a graticule shapefile that you can use as an overlay within your qgis map viewer.&lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Arial&apos;;&quot;&gt;Please enter all units in decimal degrees&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="unfinished">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Arial&apos;; font-size:11pt;&quot;&gt;&lt;span style=&quot; font-size:10pt;&quot;&gt;Tento zásuvný modul vám pomôže zostrojiť súbor shape so zemepisnou sieťou, ktorú možno preložiť cez mapy vo vašom mapovom prehliadači QGIS.&lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Arial&apos;;&quot;&gt;Všetky hodnoty zadávajte v desiatkových stupňoch&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
</context>
<context>
    <name>QgsHelpViewer</name>
    <message>
        <source>Quantum GIS Help - </source>
        <translation>Quantum GIS Pomocník -</translation>
    </message>
    <message>
        <source>Failed to get the help text from the database</source>
        <translation>Zlyhal pokus získať pomocný text z databázy</translation>
    </message>
    <message>
        <source>Error</source>
        <translation>Chyba</translation>
    </message>
    <message>
        <source>The QGIS help database is not installed</source>
        <translation>Databáza Pomocníka QGIS nie je nainštalovaná</translation>
    </message>
    <message>
        <source>This help file does not exist for your language</source>
        <translation>Tento súbor pomocníka nie je preložený do vášho jazyka</translation>
    </message>
    <message>
        <source>If you would like to create it, contact the QGIS development team</source>
        <translation>Pokiaľ aby ste radi pomohli s jeho vytvorením, kontaktujte vývojársky tím QGIS</translation>
    </message>
    <message>
        <source>Quantum GIS Help</source>
        <translation>Pomocník Quantum GIS</translation>
    </message>
</context>
<context>
    <name>QgsHelpViewerBase</name>
    <message>
        <source>QGIS Help</source>
        <translation>QGIS Pomocník</translation>
    </message>
    <message>
        <source>&amp;Home</source>
        <translation>&amp;Domov</translation>
    </message>
    <message>
        <source>Alt+H</source>
        <translation>Alt+H</translation>
    </message>
    <message>
        <source>&amp;Forward</source>
        <translation>&amp;Dopredu</translation>
    </message>
    <message>
        <source>Alt+F</source>
        <translation>Alt+F</translation>
    </message>
    <message>
        <source>&amp;Back</source>
        <translation>&amp;Späť</translation>
    </message>
    <message>
        <source>Alt+B</source>
        <translation>Alt+B</translation>
    </message>
    <message>
        <source>&amp;Close</source>
        <translation>&amp;Zatvoriť</translation>
    </message>
    <message>
        <source>Alt+C</source>
        <translation type="unfinished">Alt+Z</translation>
    </message>
</context>
<context>
    <name>QgsHttpTransaction</name>
    <message>
        <source>WMS Server responded unexpectedly with HTTP Status Code %1 (%2)</source>
        <translation>WMS Server odpovedal neočakávane s kódom stavu HTTP %1 (%2)</translation>
    </message>
    <message>
        <source>HTTP response completed, however there was an error: %1</source>
        <translation>Požiadavka HTTP dokončená, avšak vyskytla sa chyba: %1</translation>
    </message>
    <message>
        <source>Network timed out after %1 seconds of inactivity.
This may be a problem in your network connection or at the WMS server.</source>
        <translation type="obsolete">Čas požiadavky vypršal po %1 sekundách neaktivity.Problém môže byť vo vašom sieťovom spojení alebo na strane WMS servera.
        </translation>
    </message>
    <message>
        <source>HTTP transaction completed, however there was an error: %1</source>
        <translation>HTTP prenos dokončený, avšak vyskytla sa chyba: %1</translation>
    </message>
</context>
<context>
    <name>QgsIDWInterpolatorDialogBase</name>
    <message>
        <source>Dialog</source>
        <translation type="unfinished">Dialógové okno</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:12pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;span style=&quot; font-weight:600;&quot;&gt;Inverse Distance Weighting&lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-weight:600;&quot;&gt;&lt;span style=&quot; font-weight:400;&quot;&gt;The only parameter for the IDW interpolation method is the coefficient that describes the decrease of weights with distance.&lt;/span&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="unfinished">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:12pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;span style=&quot; font-weight:600;&quot;&gt;Inverse Distance Weighting&lt;/span&gt;&lt;/p&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-weight:600;&quot;&gt;&lt;span style=&quot; font-weight:400;&quot;&gt;Jediným parametrom pre interpolačnú metódu IDW je koeficient, ktorý popisuje pokles váhy so vzdialenosťou.&lt;/span&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>Distance coefficient P:</source>
        <translation type="unfinished">Koeficient vzdialenosti P:</translation>
    </message>
</context>
<context>
    <name>QgsIdentifyResults</name>
    <message>
        <source>Identify Results - </source>
        <translation>Výsledky identifikácie - </translation>
    </message>
    <message>
        <source>Run action</source>
        <translation>Spustiť akciu</translation>
    </message>
    <message>
        <source>(Derived)</source>
        <translation>(Odvodené)</translation>
    </message>
    <message>
        <source>Feature</source>
        <translation>Objekt</translation>
    </message>
    <message>
        <source>Value</source>
        <translation>Hodnota</translation>
    </message>
</context>
<context>
    <name>QgsIdentifyResultsBase</name>
    <message>
        <source>Identify Results</source>
        <translation>Výsledky identifikácie</translation>
    </message>
    <message>
        <source>Help</source>
        <translation>Pomocník</translation>
    </message>
    <message>
        <source>F1</source>
        <translation>F1</translation>
    </message>
    <message>
        <source>Close</source>
        <translation>Zatvoriť</translation>
    </message>
</context>
<context>
    <name>QgsInterpolationDialog</name>
    <message>
        <source>Triangular interpolation (TIN)</source>
        <translation type="unfinished">Triangulačná interpolácia (TIN)</translation>
    </message>
    <message>
        <source>Inverse Distance Weighting (IDW)</source>
        <translation>Metóda inverznej vzdialenosti (IDW)</translation>
    </message>
</context>
<context>
    <name>QgsInterpolationDialogBase</name>
    <message>
        <source>Interpolation plugin</source>
        <translation>Zásuvný modul na interpoláciu</translation>
    </message>
    <message>
        <source>Input</source>
        <translation>Vstup</translation>
    </message>
    <message>
        <source>Input vector layer</source>
        <translation>Vstupná vektorová vrstva</translation>
    </message>
    <message>
        <source>Use z-Coordinate for interpolation</source>
        <translation>Pre interpoláciu použiť z-ovú súradnicu </translation>
    </message>
    <message>
        <source>Interpolation attribute </source>
        <translation type="unfinished">Atribút vstupujúci do interpolácie</translation>
    </message>
    <message>
        <source>Output</source>
        <translation>Výstup</translation>
    </message>
    <message>
        <source>Interpolation method</source>
        <translation>Interpolačná metóda</translation>
    </message>
    <message>
        <source>...</source>
        <translation>...</translation>
    </message>
    <message>
        <source>Number of columns</source>
        <translation>Počet stĺpcov</translation>
    </message>
    <message>
        <source>Number of rows</source>
        <translation>Počet riadkov</translation>
    </message>
    <message>
        <source>Output file </source>
        <translation>Výstupný súbor</translation>
    </message>
</context>
<context>
    <name>QgsInterpolationPlugin</name>
    <message>
        <source>&amp;Interpolation</source>
        <translation>&amp;Interpolácia</translation>
    </message>
</context>
<context>
    <name>QgsLUDialogBase</name>
    <message>
        <source>Enter class bounds</source>
        <translation>Vložiť hranice triedy
</translation>
    </message>
    <message>
        <source>Lower value</source>
        <translation>Spodná hodnota</translation>
    </message>
    <message>
        <source>-</source>
        <translation>-</translation>
    </message>
    <message>
        <source>OK</source>
        <translation type="obsolete">OK</translation>
    </message>
    <message>
        <source>Cancel</source>
        <translation type="obsolete">Zrušiť</translation>
    </message>
    <message>
        <source>Upper value</source>
        <translation>Vrchná hodnota</translation>
    </message>
</context>
<context>
    <name>QgsLabelDialog</name>
    <message>
        <source>Auto</source>
        <translation type="unfinished">Automaticky</translation>
    </message>
</context>
<context>
    <name>QgsLabelDialogBase</name>
    <message>
        <source>Form1</source>
        <translation>Popisy</translation>
    </message>
    <message>
        <source>Field containing label:</source>
        <translation type="obsolete">Pole obsahujúce popisy:</translation>
    </message>
    <message>
        <source>Default label:</source>
        <translation type="obsolete">Štandardný text popisu:</translation>
    </message>
    <message>
        <source>Preview:</source>
        <translation>Náhľad:</translation>
    </message>
    <message>
        <source>QGIS Rocks!</source>
        <translation>QGIS Rocks!</translation>
    </message>
    <message>
        <source>Font Style</source>
        <translation type="obsolete">Štýl písma</translation>
    </message>
    <message>
        <source>Font</source>
        <translation>Písmo</translation>
    </message>
    <message>
        <source>Points</source>
        <translation>bodoch</translation>
    </message>
    <message>
        <source>Map units</source>
        <translation>mapových jednotkách</translation>
    </message>
    <message>
        <source>%</source>
        <translation>%</translation>
    </message>
    <message>
        <source>Transparency:</source>
        <translation>Priehľadnosť:</translation>
    </message>
    <message>
        <source>Colour</source>
        <translation type="obsolete">Farba</translation>
    </message>
    <message>
        <source>Position</source>
        <translation>Umiestnenie</translation>
    </message>
    <message>
        <source>X Offset (pts):</source>
        <translation type="obsolete">X-ový posun (v bodoch):</translation>
    </message>
    <message>
        <source>Y Offset (pts):</source>
        <translation type="obsolete">Y-ový posun (v bodoch):</translation>
    </message>
    <message>
        <source>Buffer Labels?</source>
        <translation type="obsolete">Okolie popisu?</translation>
    </message>
    <message>
        <source>Size:</source>
        <translation>Veľkosť:</translation>
    </message>
    <message>
        <source>Size is in map units</source>
        <translation>Veľkosť je v mapových jednotkách</translation>
    </message>
    <message>
        <source>Size is in points</source>
        <translation>Veľkosť je v bodoch</translation>
    </message>
    <message>
        <source>Above</source>
        <translation>Nad</translation>
    </message>
    <message>
        <source>Over</source>
        <translation>Skrz</translation>
    </message>
    <message>
        <source>Left</source>
        <translation>Naľavo</translation>
    </message>
    <message>
        <source>Below</source>
        <translation>Pod</translation>
    </message>
    <message>
        <source>Right</source>
        <translation>Napravo</translation>
    </message>
    <message>
        <source>Above Right</source>
        <translation>Napravo nad</translation>
    </message>
    <message>
        <source>Below Right</source>
        <translation>Napravo pod</translation>
    </message>
    <message>
        <source>Above Left</source>
        <translation>Naľavo nad</translation>
    </message>
    <message>
        <source>Below Left</source>
        <translation>Naľavo pod</translation>
    </message>
    <message>
        <source>Angle (deg):</source>
        <translation type="obsolete">Uhol (v stupňoch):</translation>
    </message>
    <message>
        <source>Data Defined Style</source>
        <translation type="obsolete">Údajmi určovaný štýl</translation>
    </message>
    <message>
        <source>&amp;Font family:</source>
        <translation type="obsolete">&amp;Rodina písiem:</translation>
    </message>
    <message>
        <source>&amp;Italic:</source>
        <translation type="obsolete">&amp;Kurzíva:</translation>
    </message>
    <message>
        <source>&amp;Underline:</source>
        <translation type="obsolete">&amp;Podčiarknuté:</translation>
    </message>
    <message>
        <source>&amp;Bold:</source>
        <translation type="obsolete">&amp;Tučné:</translation>
    </message>
    <message>
        <source>&amp;Size:</source>
        <translation type="obsolete">&amp;Veľkosť:</translation>
    </message>
    <message>
        <source>X Coordinate:</source>
        <translation type="obsolete">X-ová súradnica:</translation>
    </message>
    <message>
        <source>Y Coordinate:</source>
        <translation type="obsolete">Y-ová súradnica:</translation>
    </message>
    <message>
        <source>Placement:</source>
        <translation type="obsolete">Umiestnenie:</translation>
    </message>
    <message>
        <source>Font size units</source>
        <translation>Jednotky veľkosti písma</translation>
    </message>
    <message>
        <source>Font Alignment</source>
        <translation type="obsolete">Zarovnanie písma</translation>
    </message>
    <message>
        <source>Placement</source>
        <translation>Umiestnenie</translation>
    </message>
    <message>
        <source>Buffer</source>
        <translation>Okolie (buffer)</translation>
    </message>
    <message>
        <source>Buffer size units</source>
        <translation>Jednotky veľkosti okolia</translation>
    </message>
    <message>
        <source>Offset units</source>
        <translation>Jednotky posunutia</translation>
    </message>
    <message>
        <source>Data Defined Alignment</source>
        <translation type="obsolete">Údajmi určované zarovnanie</translation>
    </message>
    <message>
        <source>Data Defined Buffer</source>
        <translation type="obsolete">Údajmi určované okolie</translation>
    </message>
    <message>
        <source>Data Defined Position</source>
        <translation type="obsolete">Údajmi určovaná poloha</translation>
    </message>
    <message>
        <source>Source</source>
        <translation type="obsolete">Zdroj</translation>
    </message>
    <message>
        <source>Size Units:</source>
        <translation type="obsolete">Jednotky veľkosti:</translation>
    </message>
    <message>
        <source>Field containing label</source>
        <translation type="unfinished">Pole obsahujúce popis</translation>
    </message>
    <message>
        <source>Default label</source>
        <translation type="unfinished">Štandardný text popisu</translation>
    </message>
    <message>
        <source>Data defined style</source>
        <translation type="unfinished">Údajmi určovaný štýl</translation>
    </message>
    <message>
        <source>Data defined alignment</source>
        <translation type="unfinished">Údajmi určované zarovnanie</translation>
    </message>
    <message>
        <source>Data defined buffer</source>
        <translation type="unfinished">Údajmi určované okolie</translation>
    </message>
    <message>
        <source>Data defined position</source>
        <translation type="unfinished">Údajmi určovaná poloha</translation>
    </message>
    <message>
        <source>Font transparency</source>
        <translation type="unfinished">Priehľadnosť písma</translation>
    </message>
    <message>
        <source>Color</source>
        <translation>Farba</translation>
    </message>
    <message>
        <source>Angle (deg)</source>
        <translation>Uhol (v stupňoch)</translation>
    </message>
    <message>
        <source>&#xb0;</source>
        <translation type="obsolete">°</translation>
    </message>
    <message>
        <source>Buffer labels?</source>
        <translation>Okolie popisu?</translation>
    </message>
    <message>
        <source>Buffer size</source>
        <translation>Veľkosť okolia</translation>
    </message>
    <message>
        <source>Transparency</source>
        <translation>Priehľadnosť</translation>
    </message>
    <message>
        <source>Multiline labels?</source>
        <translation>Viacriadkové popisy?</translation>
    </message>
    <message>
        <source>X Offset (pts)</source>
        <translation>X-ový posun (v bodoch)</translation>
    </message>
    <message>
        <source>Y Offset (pts)</source>
        <translation>Y-ový posun (v bodoch)</translation>
    </message>
    <message>
        <source>&amp;Font family</source>
        <translation>&amp;Rodina písiem</translation>
    </message>
    <message>
        <source>&amp;Bold</source>
        <translation>&amp;Tučné</translation>
    </message>
    <message>
        <source>&amp;Italic</source>
        <translation>&amp;Kurzíva</translation>
    </message>
    <message>
        <source>&amp;Underline</source>
        <translation>&amp;Podčiarknuté</translation>
    </message>
    <message>
        <source>&amp;Size</source>
        <translation>&amp;Veľkosť</translation>
    </message>
    <message>
        <source>Size units</source>
        <translation>Jednotky veľkosti</translation>
    </message>
    <message>
        <source>X Coordinate</source>
        <translation>X-ová súradnica</translation>
    </message>
    <message>
        <source>Y Coordinate</source>
        <translation>Y-ová súradnica</translation>
    </message>
    <message>
        <source>General</source>
        <translation>Všeobecné</translation>
    </message>
    <message>
        <source>Use scale dependent rendering</source>
        <translation>Používať vykresľovanie v závislosti od mierky</translation>
    </message>
    <message>
        <source>Maximum</source>
        <translation>Maximum</translation>
    </message>
    <message>
        <source>Minimum</source>
        <translation>Minimum</translation>
    </message>
    <message>
        <source>Minimum scale at which this layer will be displayed. </source>
        <translation>Minimálna mierka pri ktorej bude táto vrstva zobrazená.</translation>
    </message>
    <message>
        <source>Maximum scale at which this layer will be displayed. </source>
        <translation>Maximálna mierka pri ktorej bude táto vrstva zobrazená.</translation>
    </message>
</context>
<context>
    <name>QgsLayerProjectionSelector</name>
    <message>
        <source>Define this layer&apos;s projection:</source>
        <translation type="obsolete">Určiť mapové zobrazenie tejto vrstvy:</translation>
    </message>
    <message>
        <source>This layer appears to have no projection specification.</source>
        <translation type="obsolete">Táto vrstva zrejme nemá uvedenú žiadnu informáciu o použitom mapovom zobrazení.</translation>
    </message>
    <message>
        <source>By default, this layer will now have its projection set to that of the project, but you may override this by selecting a different projection below.</source>
        <translation type="obsolete">Predvolene bude mať teraz nastavené mapové zobrazenie zhodné s nastavením projektu, ale je možné ho zmeniť vybratím iného zobraznia z nižšie uvedenej ponuky.</translation>
    </message>
</context>
<context>
    <name>QgsLayerProjectionSelectorBase</name>
    <message>
        <source>Layer Projection Selector</source>
        <translation type="obsolete">Výber mapového zobrazenia vrstvy</translation>
    </message>
    <message>
        <source>OK</source>
        <translation type="obsolete">OK</translation>
    </message>
</context>
<context>
    <name>QgsLegend</name>
    <message>
        <source>group</source>
        <translation>skupina</translation>
    </message>
    <message>
        <source>&amp;Remove</source>
        <translation>&amp;Odobrať</translation>
    </message>
    <message>
        <source>&amp;Make to toplevel item</source>
        <translation>&amp;Premiestniť položku do najvyššej úrovne</translation>
    </message>
    <message>
        <source>Re&amp;name</source>
        <translation>Preme&amp;novať</translation>
    </message>
    <message>
        <source>&amp;Add group</source>
        <translation>Prid&amp;ať skupinu</translation>
    </message>
    <message>
        <source>&amp;Expand all</source>
        <translation>&amp;Rozbaliť</translation>
    </message>
    <message>
        <source>&amp;Collapse all</source>
        <translation>&amp;Zabaliť</translation>
    </message>
    <message>
        <source>Show file groups</source>
        <translation>Ukázať skupiny súborov</translation>
    </message>
    <message>
        <source>No Layer Selected</source>
        <translation type="unfinished">Nie je vybratá žiadna vrstva</translation>
    </message>
    <message>
        <source>To open an attribute table, you must select a vector layer in the legend</source>
        <translation type="unfinished">Pred otvorením tabuľky atribútov je nutné vybrať vrstvu v okne Legenda</translation>
    </message>
</context>
<context>
    <name>QgsLegendLayer</name>
    <message>
        <source>&amp;Zoom to layer extent</source>
        <translation type="unfinished">Pohľad na veľkosť &amp;vrstvy</translation>
    </message>
    <message>
        <source>&amp;Zoom to best scale (100%)</source>
        <translation type="unfinished">&amp;Zmeniť pohľad na 1:1 (100%)</translation>
    </message>
    <message>
        <source>&amp;Show in overview</source>
        <translation type="unfinished">&amp;Ukázať v prehľade</translation>
    </message>
    <message>
        <source>&amp;Remove</source>
        <translation>&amp;Odobrať</translation>
    </message>
    <message>
        <source>&amp;Open attribute table</source>
        <translation type="unfinished">Otvoriť &amp;tabuľku atribútov</translation>
    </message>
    <message>
        <source>Save as shapefile...</source>
        <translation>Uložiť ako shape súbor...</translation>
    </message>
    <message>
        <source>Save selection as shapefile...</source>
        <translation>Uložiť výber ako súbor shape...</translation>
    </message>
    <message>
        <source>&amp;Properties</source>
        <translation>&amp;Vlastnosti</translation>
    </message>
    <message>
        <source>More layers</source>
        <translation type="obsolete">Viac vrstiev</translation>
    </message>
    <message>
        <source>This item contains more layer files. Displaying more layers in table is not supported.</source>
        <translation type="obsolete">Táto položka obsahuje viac súborov s vrstvami. Zobrazenie viacerých vrstiev v tabuľke nie je podporované.</translation>
    </message>
    <message>
        <source>Multiple layers</source>
        <translation>Viacero vrstiev</translation>
    </message>
    <message>
        <source>This item contains multiple layers. Displaying multiple layers in the table is not supported.</source>
        <translation>Táto položka obsahuje viacero vrstiev. Zobrazenie viacerých vrstiev v tabuľke nie je podporované.</translation>
    </message>
</context>
<context>
    <name>QgsLegendLayerFile</name>
    <message>
        <source>Not a vector layer</source>
        <translation type="obsolete">Nie je vektorová vrstva</translation>
    </message>
    <message>
        <source>To open an attribute table, you must select a vector layer in the legend</source>
        <translation type="obsolete">Pred otvorením tabuľky atribútov je nutné vybrať vrstvu v okne Legenda</translation>
    </message>
    <message>
        <source>Attribute table - </source>
        <translation type="obsolete">Tabuľka atribútov -</translation>
    </message>
    <message>
        <source>Save layer as...</source>
        <translation>Uložiť vrstvu ako...</translation>
    </message>
    <message>
        <source>Saving done</source>
        <translation type="unfinished">Ukladanie dokončené</translation>
    </message>
    <message>
        <source>Export to Shapefile has been completed</source>
        <translation type="unfinished">Export do súboru shape bol dokončený</translation>
    </message>
    <message>
        <source>Driver not found</source>
        <translation>Ovládač sa nenašiel</translation>
    </message>
    <message>
        <source>ESRI Shapefile driver is not available</source>
        <translation type="unfinished">Ovládač pre formát ESRI shape nie je dostupný</translation>
    </message>
    <message>
        <source>Error creating shapefile</source>
        <translation type="unfinished">Chyba pri vytváraní súboru Shape</translation>
    </message>
    <message>
        <source>The shapefile could not be created (</source>
        <translation type="unfinished">Súbor Shape nemožno vytvoriť (</translation>
    </message>
    <message>
        <source>Error</source>
        <translation>Chyba</translation>
    </message>
    <message>
        <source>Layer creation failed</source>
        <translation type="unfinished">Vytvorenie vrstvy zlyhalo</translation>
    </message>
    <message>
        <source>Start editing failed</source>
        <translation type="obsolete">Pokus o úpravy zlyhal</translation>
    </message>
    <message>
        <source>Provider cannot be opened for editing</source>
        <translation type="obsolete">Správca sa nedá otvoriť pre zápis</translation>
    </message>
    <message>
        <source>Stop editing</source>
        <translation type="obsolete">Ukončiť úpravy</translation>
    </message>
    <message>
        <source>Do you want to save the changes?</source>
        <translation type="obsolete">Prajete si uložiť zmeny?</translation>
    </message>
    <message>
        <source>Could not commit changes</source>
        <translation type="obsolete">Nemožno odoslať zmeny </translation>
    </message>
    <message>
        <source>Problems during roll back</source>
        <translation type="obsolete">Problémy v priebehu návratu do východzieho stavu (roll back)</translation>
    </message>
    <message>
        <source>&amp;Zoom to layer extent</source>
        <translation type="unfinished">Pohľad na veľkosť &amp;vrstvy</translation>
    </message>
    <message>
        <source>&amp;Show in overview</source>
        <translation type="unfinished">&amp;Ukázať v prehľade</translation>
    </message>
    <message>
        <source>&amp;Remove</source>
        <translation>&amp;Odobrať</translation>
    </message>
    <message>
        <source>&amp;Open attribute table</source>
        <translation>Otvoriť &amp;tabuľku atribútov</translation>
    </message>
    <message>
        <source>Save as shapefile...</source>
        <translation>Uložiť ako shape súbor...</translation>
    </message>
    <message>
        <source>Save selection as shapefile...</source>
        <translation>Uložiť výber ako súbor shape...</translation>
    </message>
    <message>
        <source>&amp;Properties</source>
        <translation>&amp;Vlastnosti</translation>
    </message>
    <message>
        <source>bad_alloc exception</source>
        <translation type="obsolete">výnimka bad_alloc</translation>
    </message>
    <message>
        <source>Filling the attribute table has been stopped because there was no more virtual memory left</source>
        <translation type="obsolete">Naplnenie tabuľky atribútov bolo zastavené, pretože už nebol dostatok virtuálnej pamäte</translation>
    </message>
    <message>
        <source>Layer attribute table contains unsupported datatype(s)</source>
        <translation>Atribútová tabuľka vrstvy obsahuje napodporovaný(é) typ(y) údajov</translation>
    </message>
    <message>
        <source>Select the coordinate reference system for the saved shapefile.</source>
        <translation>Vyberte referenčný súradnicový systém pre ukladaný súbor vo formáte shape.</translation>
    </message>
    <message>
        <source>The data points will be transformed from the layer coordinate reference system.</source>
        <translation>Údajové body budú prevedené z referenčného súradnicového systému vrstvy.</translation>
    </message>
</context>
<context>
    <name>QgsLineStyleDialogBase</name>
    <message>
        <source>Select a line style</source>
        <translation type="obsolete">Výber štýlu pre línie</translation>
    </message>
    <message>
        <source>Styles</source>
        <translation type="obsolete">Štýly</translation>
    </message>
    <message>
        <source>Ok</source>
        <translation type="obsolete">OK</translation>
    </message>
    <message>
        <source>Cancel</source>
        <translation type="obsolete">Zrušiť</translation>
    </message>
</context>
<context>
    <name>QgsLineStyleWidgetBase</name>
    <message>
        <source>Form2</source>
        <translation type="obsolete">Štýl línie</translation>
    </message>
    <message>
        <source>Outline Style</source>
        <translation type="obsolete">Štýl obrysu</translation>
    </message>
    <message>
        <source>Width:</source>
        <translation type="obsolete">Hrúbka:</translation>
    </message>
    <message>
        <source>Colour:</source>
        <translation type="obsolete">Farba:</translation>
    </message>
    <message>
        <source>LineStyleWidget</source>
        <translation type="obsolete">LineStyleWidget</translation>
    </message>
    <message>
        <source>col</source>
        <translation type="obsolete">col</translation>
    </message>
</context>
<context>
    <name>QgsMapCanvas</name>
    <message>
        <source>Could not draw</source>
        <translation>Nemožno vykresľovať</translation>
    </message>
    <message>
        <source>because</source>
        <translation>z nasledovného dôvodu</translation>
    </message>
</context>
<context>
    <name>QgsMapLayer</name>
    <message>
        <source> Check file permissions and retry.</source>
        <translation type="obsolete"> Skontrolujte nastavenia prístupových práv a skúste znova.</translation>
    </message>
    <message>
        <source>%1 at line %2 column %3</source>
        <translation>%1 na riadku %2 stĺpec %3</translation>
    </message>
    <message>
        <source>style not found in database</source>
        <translation>štýl sa nenšiel v databáze</translation>
    </message>
    <message>
        <source>User database could not be opened.</source>
        <translation>Nemožno otvoriť užívateľskú databázu.</translation>
    </message>
    <message>
        <source>The style table could not be created.</source>
        <translation>Nemožno vytvoriť tabuľku so štýlom.</translation>
    </message>
    <message>
        <source>The style %1 was saved to database</source>
        <translation>Štýl %1 bol uložený do databázy</translation>
    </message>
    <message>
        <source>The style %1 was updated in the database.</source>
        <translation>Štýl %1 bol aktualizovaný v databáze.</translation>
    </message>
    <message>
        <source>The style %1 could not be updated in the database.</source>
        <translation>Štýl %1 sa v databáze nepodarilo aktualizovať.</translation>
    </message>
    <message>
        <source>The style %1 could not be inserted into database.</source>
        <translation>Štýl %1 nemožno vložiť do databázy.</translation>
    </message>
</context>
<context>
    <name>QgsMapToolIdentify</name>
    <message>
        <source>No features found</source>
        <translation type="obsolete">Nenašli sa žiadne objekty</translation>
    </message>
    <message>
        <source>&lt;p&gt;No features were found within the search radius. Note that it is currently not possible to use the identify tool on unsaved features.&lt;/p&gt;</source>
        <translation type="obsolete">&lt;p&gt;V okruhu danom poloemor vyhľadávania sa nenašli žiadne objekty.&lt;/p&gt;&lt;p&gt;Poznámka: V súčasnosti nie je možné využívať nástroj identifikácie na neuložené objekty.&lt;/p&gt;</translation>
    </message>
    <message>
        <source>- %1 features found</source>
        <comment>Identify results window title</comment>
        <translation type="obsolete">- nájdených %1 objektov
        </translation>
    </message>
    <message>
        <source>(clicked coordinate)</source>
        <translation>(kliknutá súradnica)</translation>
    </message>
    <message>
        <source>WMS identify result for %1
%2</source>
        <translation>Výsledok identifikácie WMS pre %1
%2</translation>
    </message>
</context>
<context>
    <name>QgsMapToolSplitFeatures</name>
    <message>
        <source>Split error</source>
        <translation>Chyba pri rozdeľovaní</translation>
    </message>
    <message>
        <source>An error occured during feature splitting</source>
        <translation>Pri rozdeľovaní objektu nastala chyba</translation>
    </message>
    <message>
        <source>No feature split done</source>
        <translation>Objekt nebol rozdelený</translation>
    </message>
    <message>
        <source>If there are selected features, the split tool only applies to the selected ones. If you like to split all features under the split line, clear the selection</source>
        <translation>Ak sú vybrané objekty, nástroj na rozdelenie bude použitý len na tieto vybrané objekty. Pokiaľ majú byť rozdelené všetky objekty pod líniou rozdelenia, je potrebné vyčistiť výber</translation>
    </message>
</context>
<context>
    <name>QgsMapToolVertexEdit</name>
    <message>
        <source>Snap tolerance</source>
        <translation>Tolerancia zameriavania</translation>
    </message>
    <message>
        <source>Don&apos;t show this message again</source>
        <translation>Nezobrazovať nabudúce túto správu</translation>
    </message>
    <message>
        <source>Could not snap segment.</source>
        <translation>Nemožno zamerať na úsek.</translation>
    </message>
    <message>
        <source>Have you set the tolerance in Settings &gt; Project Properties &gt; General?</source>
        <translation>Je správne nastavená tolerancia v menu Nastavenia &gt; Vlastnosti projektu &gt; Všeobecné?</translation>
    </message>
</context>
<context>
    <name>QgsMapserverExport</name>
    <message>
        <source>Overwrite File?</source>
        <translation>Prepísať súbor?</translation>
    </message>
    <message>
        <source> exists. 
Do you want to overwrite it?</source>
        <translation> existuje. 
Želáte si ho prepísať?</translation>
    </message>
    <message>
        <source>Name for the map file</source>
        <translation type="unfinished">Meno pre súbor map</translation>
    </message>
    <message>
        <source>Choose the QGIS project file</source>
        <translation>Vyberte súbor QGIS projektu</translation>
    </message>
    <message>
        <source>QGIS Project Files (*.qgs);;All files (*.*)</source>
        <comment>Filter list for selecting files from a dialog box</comment>
        <translation>Súbory projektov QGIS (*.qgs);;Všetky súbory (*.*)</translation>
    </message>
    <message>
        <source> exists. 
Do you want to overwrite it?</source>
        <comment>a filename is prepended to this text, and appears in a dialog box</comment>
        <translation type="obsolete"> existuje. Želáte si ho prepísať?</translation>
    </message>
    <message>
        <source>MapServer map files (*.map);;All files (*.*)</source>
        <comment>

Filter list for selecting files from a dialog box</comment>
        <translation type="obsolete">Súbory map pre MapServer (*.map);;Všetky súbory(*.*)</translation>
    </message>
    <message>
        <source>MapServer map files (*.map);;All files (*.*)</source>
        <comment>Filter list for selecting files from a dialog box</comment>
        <translation>Súbory map pre MapServer (*.map);;Všetky súbory(*.*)</translation>
    </message>
    <message>
        <source> exists. 
Do you want to overwrite it?</source>
        <comment>a fileName is prepended to this text, and appears in a dialog box</comment>
        <translation> už existuje. 
Želáte si ho prepísať?</translation>
    </message>
</context>
<context>
    <name>QgsMapserverExportBase</name>
    <message>
        <source>Export to Mapserver</source>
        <translation>Exportovať ako Mapserver</translation>
    </message>
    <message>
        <source>Map file</source>
        <translation>Mapový súbor</translation>
    </message>
    <message>
        <source>Export LAYER information only</source>
        <translation>Exportovať len informáciu o VRSTVE</translation>
    </message>
    <message>
        <source>Map</source>
        <translation>Mapa</translation>
    </message>
    <message>
        <source>Name</source>
        <translation>Meno</translation>
    </message>
    <message>
        <source>Height</source>
        <translation>Výška</translation>
    </message>
    <message>
        <source>Width</source>
        <translation>Šírka</translation>
    </message>
    <message>
        <source>dd</source>
        <translation>dd</translation>
    </message>
    <message>
        <source>feet</source>
        <translation>feet</translation>
    </message>
    <message>
        <source>meters</source>
        <translation>meters</translation>
    </message>
    <message>
        <source>miles</source>
        <translation>miles</translation>
    </message>
    <message>
        <source>inches</source>
        <translation>inches</translation>
    </message>
    <message>
        <source>kilometers</source>
        <translation>kilometers</translation>
    </message>
    <message>
        <source>Units</source>
        <translation>Jednotky</translation>
    </message>
    <message>
        <source>Image type</source>
        <translation>Typ obrázku</translation>
    </message>
    <message>
        <source>gif</source>
        <translation>gif</translation>
    </message>
    <message>
        <source>gtiff</source>
        <translation>gtiff</translation>
    </message>
    <message>
        <source>jpeg</source>
        <translation>jpeg</translation>
    </message>
    <message>
        <source>png</source>
        <translation>png</translation>
    </message>
    <message>
        <source>swf</source>
        <translation>swf</translation>
    </message>
    <message>
        <source>userdefined</source>
        <translation>userdefined</translation>
    </message>
    <message>
        <source>wbmp</source>
        <translation>wbmp</translation>
    </message>
    <message>
        <source>MinScale</source>
        <translation>MinMierka</translation>
    </message>
    <message>
        <source>MaxScale</source>
        <translation>MaxMierka</translation>
    </message>
    <message>
        <source>Prefix attached to map, scalebar and legend GIF filenames created using this MapFile. It should be kept short.</source>
        <translation>Predpona pripojená k názvom GIF súborov mapy, grafickej mierky a legendy vytvorených použitím tohoto mapového súboru. Mala by byť krátka.</translation>
    </message>
    <message>
        <source>Web Interface Definition</source>
        <translation>Definícia webového rozhrania</translation>
    </message>
    <message>
        <source>Header</source>
        <translation>Hlavička</translation>
    </message>
    <message>
        <source>Footer</source>
        <translation>Pätička</translation>
    </message>
    <message>
        <source>Template</source>
        <translation>Šablóna</translation>
    </message>
    <message>
        <source>&amp;Help</source>
        <translation>&amp;Pomocník</translation>
    </message>
    <message>
        <source>F1</source>
        <translation>F1</translation>
    </message>
    <message>
        <source>&amp;OK</source>
        <translation>&amp;OK</translation>
    </message>
    <message>
        <source>&amp;Cancel</source>
        <translation>&amp;Zrušiť</translation>
    </message>
    <message>
        <source>...</source>
        <translation>...</translation>
    </message>
    <message>
        <source>Name for the map file to be created from the QGIS project file</source>
        <translation>Meno mapového súboru, ktorý má byť vytvorený na základe súboru projektu pre QGIS </translation>
    </message>
    <message>
        <source>If checked, only the layer information will be processed</source>
        <translation>Pokiaľ je zaškrtnuté, bude len spracovaná informácia o vrstve</translation>
    </message>
    <message>
        <source>Path to the MapServer template file</source>
        <translation>Cesta k súboru so šablónou pre MapServer</translation>
    </message>
    <message>
        <source>Prefix attached to map, scalebar and legend GIF filenames created using this MapFile</source>
        <translation type="unfinished">Predpona pripojená k mape, mená súborov GIS s grafickou mierkou a legendou vytvorené s použitím tohoto MapFile</translation>
    </message>
    <message>
        <source>Full path to the QGIS project file to export to MapServer map format</source>
        <translation>Úplná cesta k súboru projektu QGIS určeného na export do formátu map pre MapServer</translation>
    </message>
    <message>
        <source>QGIS project file</source>
        <translation>Súbor vo formáte projektu QGIS</translation>
    </message>
    <message>
        <source>Browse...</source>
        <translation>Prechádzať...</translation>
    </message>
    <message>
        <source>Save As...</source>
        <translation>Uložiť ako...</translation>
    </message>
</context>
<context>
    <name>QgsMarkerDialogBase</name>
    <message>
        <source>Choose a marker symbol</source>
        <translation type="obsolete">Vybrať symbol pre značku</translation>
    </message>
    <message>
        <source>Directory</source>
        <translation type="obsolete">Adresár</translation>
    </message>
    <message>
        <source>...</source>
        <translation type="obsolete">...</translation>
    </message>
    <message>
        <source>Ok</source>
        <translation type="obsolete">OK</translation>
    </message>
    <message>
        <source>Cancel</source>
        <translation type="obsolete">Zrušiť</translation>
    </message>
    <message>
        <source>New Item</source>
        <translation type="obsolete">Nová položka</translation>
    </message>
</context>
<context>
    <name>QgsMeasureBase</name>
    <message>
        <source>Measure</source>
        <translation>Vzdialenosť</translation>
    </message>
    <message>
        <source>New</source>
        <translation>Nové</translation>
    </message>
    <message>
        <source>Help</source>
        <translation>Pomocník</translation>
    </message>
    <message>
        <source>Cl&amp;ose</source>
        <translation>&amp;Zatvoriť</translation>
    </message>
    <message>
        <source>Total:</source>
        <translation>Spolu:</translation>
    </message>
    <message>
        <source>Segments</source>
        <translation>Úseky</translation>
    </message>
</context>
<context>
    <name>QgsMeasureDialog</name>
    <message>
        <source>Segments (in meters)</source>
        <translation>Úseky (v metroch)</translation>
    </message>
    <message>
        <source>Segments (in feet)</source>
        <translation>Úseky (v stopách)</translation>
    </message>
    <message>
        <source>Segments (in degrees)</source>
        <translation>Úseky (v stupňoch)</translation>
    </message>
    <message>
        <source>Segments</source>
        <translation>Úseky</translation>
    </message>
</context>
<context>
    <name>QgsMeasureTool</name>
    <message>
        <source>Incorrect measure results</source>
        <translation>Nesprávny výsledok merania</translation>
    </message>
    <message>
        <source>&lt;p&gt;This map is defined with a geographic coordinate system (latitude/longitude) but the map extents suggests that it is actually a projected coordinate system (e.g., Mercator). If so, the results from line or area measurements will be incorrect.&lt;/p&gt;&lt;p&gt;To fix this, explicitly set an appropriate map coordinate system using the &lt;tt&gt;Settings:Project Properties&lt;/tt&gt; menu.</source>
        <translation>Táto mapa je definovaná so zemepisným súradnicovým systémom (šírka/dĺžka) ale z rozsahu mapy vyplýva, že využíva mapový súradnicový systém (napr. Mercator). Ak je to naozaj tak, výsledky z merania dĺžky a rozlohy nebudú správne.&lt;/p&gt;&lt;p&gt;To možno napraviť zadaním správneho súradnicového systému cez menu &lt;tt&gt;Nastavenia: Nastavenia projektu&lt;/tt&gt;.</translation>
    </message>
</context>
<context>
    <name>QgsMessageViewer</name>
    <message>
        <source>QGIS Message</source>
        <translation>Správa QGIS</translation>
    </message>
    <message>
        <source>Close</source>
        <translation>Zatvoriť</translation>
    </message>
    <message>
        <source>Don&apos;t show this message again</source>
        <translation>Nezozbrazovať viac toto upozornenie</translation>
    </message>
</context>
<context>
    <name>QgsMySQLProvider</name>
    <message>
        <source>Unable to access relation</source>
        <translation type="obsolete">Nie je možné pristúpiť k relácii</translation>
    </message>
    <message>
        <source>Unable to access the </source>
        <translation type="obsolete"> Nemožno pristúpiť k relácii </translation>
    </message>
    <message>
        <source> relation.
The error message from the database was:
</source>
        <translation type="obsolete"> .
Chybové hlásenie z databázy:
</translation>
    </message>
    <message>
        <source>No GEOS Support!</source>
        <translation type="obsolete">Bez podpory GEOS!</translation>
    </message>
    <message>
        <source>Your PostGIS installation has no GEOS support.
Feature selection and identification will not work properly.
Please install PostGIS with GEOS support (http://geos.refractions.net)</source>
        <translation type="obsolete">Vaša inštalácia PostGIS nemá podporu GEOSu.
Výber objektov a identifikácia nebudú pracovať správne.
Prosím nainštalujte PostGIS s podporou GEOSu (http://geos.refractions.net)</translation>
    </message>
</context>
<context>
    <name>QgsNewConnection</name>
    <message>
        <source>Test connection</source>
        <translation>Vyskúšať spojenie</translation>
    </message>
    <message>
        <source>Connection failed - Check settings and try again.

Extended error information:
</source>
        <translation>Spojenie zlyhalo - skontrolujte nastavenia a skúste znova.

Rozšírené informácie o chybe:
</translation>
    </message>
    <message>
        <source>Connection to %1 was successful</source>
        <translation>Spojenie k databáze %1 bolo úspešné</translation>
    </message>
</context>
<context>
    <name>QgsNewConnectionBase</name>
    <message>
        <source>Create a New PostGIS connection</source>
        <translation>Vytvoriť nové spojenie PostGIS</translation>
    </message>
    <message>
        <source>OK</source>
        <translation>OK</translation>
    </message>
    <message>
        <source>Cancel</source>
        <translation>Zrušiť</translation>
    </message>
    <message>
        <source>Help</source>
        <translation>Pomocník</translation>
    </message>
    <message>
        <source>Connection Information</source>
        <translation>Informácie o spojení</translation>
    </message>
    <message>
        <source>Host</source>
        <translation>Hostiteľ</translation>
    </message>
    <message>
        <source>Database</source>
        <translation>Databáza</translation>
    </message>
    <message>
        <source>Username</source>
        <translation>Meno používateľa</translation>
    </message>
    <message>
        <source>Name</source>
        <translation>Meno</translation>
    </message>
    <message>
        <source>Name of the new connection</source>
        <translation>Názov nového spojenia</translation>
    </message>
    <message>
        <source>Password</source>
        <translation>Heslo</translation>
    </message>
    <message>
        <source>Test Connect</source>
        <translation>Vyskúšať spojenie</translation>
    </message>
    <message>
        <source>Save Password</source>
        <translation>Uložiť heslo</translation>
    </message>
    <message>
        <source>F1</source>
        <translation>F1</translation>
    </message>
    <message>
        <source>Port</source>
        <translation>Port</translation>
    </message>
    <message>
        <source>5432</source>
        <translation>5432</translation>
    </message>
    <message>
        <source>Only look in the &apos;public&apos; schema</source>
        <translation>Prezerať len v schéme &apos;public&apos;</translation>
    </message>
    <message>
        <source>Only look in the geometry_columns table</source>
        <translation>Prezerať len tabuľku geometry_columns</translation>
    </message>
    <message>
        <source>Restrict the search to the public schema for spatial tables not in the geometry_columns table</source>
        <translation type="unfinished">Pri vyhľadávaní priestorových tabuliek, nenachádzajúcich sa v tabuľke geometry_columns obmedziť vyhľadávanie len na schému public</translation>
    </message>
    <message>
        <source>When searching for spatial tables that are not in the geometry_columns tables, restrict the search to tables that are in the public schema (for some databases this can save lots of time)</source>
        <translation type="unfinished">Pri prehľadávaní priestorových tabuliek, ktoré nie sú v tabuľke geometry_columns bude hľadanie zúžené len na tabuľky ktoré sa nachádzajú v schéme public (pri neiktorých databázach to môže ušetriť veľa času)</translation>
    </message>
    <message>
        <source>Restrict the displayed tables to those that are in the geometry_columns table</source>
        <translation>Obmedziť zobrazenie tabuliek len na tie, ktoré sú uvedené v tabuľke geometry_columns</translation>
    </message>
    <message>
        <source>Restricts the displayed tables to those that are in the geometry_columns table. This can speed up the initial display of spatial tables.</source>
        <translation type="unfinished">Zúži zobrazenie tabuleik len na tie, ktoré sú uvedené v tabuľke geometry_columns. To môže urýchliť prvotné zobrazenie priestorových tabuliek.</translation>
    </message>
</context>
<context>
    <name>QgsNewHttpConnectionBase</name>
    <message>
        <source>OK</source>
        <translation type="obsolete">OK</translation>
    </message>
    <message>
        <source>Cancel</source>
        <translation type="obsolete">Zrušiť</translation>
    </message>
    <message>
        <source>Help</source>
        <translation type="obsolete">Pomocník</translation>
    </message>
    <message>
        <source>F1</source>
        <translation type="obsolete">F1</translation>
    </message>
    <message>
        <source>Connection Information</source>
        <translation type="obsolete">Informácie o spojení</translation>
    </message>
    <message>
        <source>Name</source>
        <translation>Meno</translation>
    </message>
    <message>
        <source>URL</source>
        <translation>URL</translation>
    </message>
    <message>
        <source>Name of the new connection</source>
        <translation>Meno nového spojenia</translation>
    </message>
    <message>
        <source>Create a New WMS connection</source>
        <translation type="obsolete">Vytvoriť nové spojenie k WMS</translation>
    </message>
    <message>
        <source>Proxy Host</source>
        <translation type="obsolete">Proxy server</translation>
    </message>
    <message>
        <source>Proxy Port</source>
        <translation type="obsolete">Proxy port</translation>
    </message>
    <message>
        <source>HTTP address of the Web Map Server</source>
        <translation>HTTP adresa WMS servera</translation>
    </message>
    <message>
        <source>Name of your HTTP proxy (optional)</source>
        <translation type="obsolete">Názov vášho HTTP proxy servera (voliteľné)</translation>
    </message>
    <message>
        <source>Port number of your HTTP proxy (optional)</source>
        <translation type="obsolete">Číslo portu vášho HTTP proxy servera (voliteľné)</translation>
    </message>
    <message>
        <source>Proxy User</source>
        <translation type="obsolete">Užívateľ proxy</translation>
    </message>
    <message>
        <source>Proxy Password</source>
        <translation type="obsolete">Heslo proxy</translation>
    </message>
    <message>
        <source>Your user name for the HTTP proxy (optional)</source>
        <translation type="obsolete">Vaše prihlasovacie meno pre HTTP proxy server (voliteľné)</translation>
    </message>
    <message>
        <source>Password for your HTTP proxy (optional)</source>
        <translation type="obsolete">Heslo pre váš HTTP proxy server (voliteľné)</translation>
    </message>
    <message>
        <source>Create a new WMS connection</source>
        <translation type="unfinished">Vytvoriť nové spojenie k WMS</translation>
    </message>
    <message>
        <source>Connection details</source>
        <translation type="unfinished">Detaily spojenia</translation>
    </message>
</context>
<context>
    <name>QgsNorthArrowPlugin</name>
    <message>
        <source>Bottom Left</source>
        <translation>Vľavo dole</translation>
    </message>
    <message>
        <source>Top Right</source>
        <translation>Vpravo hore</translation>
    </message>
    <message>
        <source>Bottom Right</source>
        <translation>Vpravo dole</translation>
    </message>
    <message>
        <source>&amp;Decorations</source>
        <translation>&amp;Doplnky</translation>
    </message>
    <message>
        <source>Creates a north arrow that is displayed on the map canvas</source>
        <translation>Vytvorí na mapovom plátne smerovú ružicu</translation>
    </message>
    <message>
        <source>Top Left</source>
        <translation>Vľavo hore</translation>
    </message>
    <message>
        <source>&amp;North Arrow</source>
        <translation>&amp;Smerová ružica</translation>
    </message>
    <message>
        <source>North arrow pixmap not found</source>
        <translation>Obrázok smerovej ružice sa nenašiel</translation>
    </message>
</context>
<context>
    <name>QgsNorthArrowPluginGui</name>
    <message>
        <source>Pixmap not found</source>
        <translation>Obrázok sa nenašiel</translation>
    </message>
</context>
<context>
    <name>QgsNorthArrowPluginGuiBase</name>
    <message>
        <source>North Arrow Plugin</source>
        <translation>Zásuvný modul Smerová ružica</translation>
    </message>
    <message>
        <source>Properties</source>
        <translation>Vlastnosti</translation>
    </message>
    <message>
        <source>Angle</source>
        <translation>Uhol</translation>
    </message>
    <message>
        <source>Enable North Arrow</source>
        <translation>Zapnúť smerovú ružicu</translation>
    </message>
    <message>
        <source>Top Left</source>
        <translation>Vľavo hore</translation>
    </message>
    <message>
        <source>Top Right</source>
        <translation>Vpravo hore</translation>
    </message>
    <message>
        <source>Bottom Left</source>
        <translation>Vľavo dole</translation>
    </message>
    <message>
        <source>Bottom Right</source>
        <translation>Vpravo dole</translation>
    </message>
    <message>
        <source>Placement on screen</source>
        <translation>Umiestnenie na obrazovke</translation>
    </message>
    <message>
        <source>Preview of north arrow</source>
        <translation>Náhľad na smerovú ružicu</translation>
    </message>
    <message>
        <source>Placement</source>
        <translation>Umiestnenie</translation>
    </message>
    <message>
        <source>Icon</source>
        <translation>Ikona</translation>
    </message>
    <message>
        <source>Set direction automatically</source>
        <translation>Nastaviť smer automaticky</translation>
    </message>
    <message>
        <source>Browse...</source>
        <translation>Prechádzať...</translation>
    </message>
</context>
<context>
    <name>QgsOptions</name>
    <message>
        <source>Detected active locale on your system: </source>
        <translation type="unfinished">Nájdené aktívne regionálne nastavenie na vašom systéme: </translation>
    </message>
    <message>
        <source>to vertex</source>
        <translation>k uzlu</translation>
    </message>
    <message>
        <source>to segment</source>
        <translation>k úseku</translation>
    </message>
    <message>
        <source>to vertex and segment</source>
        <translation>k uzlu a úseku</translation>
    </message>
    <message>
        <source>Semi transparent circle</source>
        <translation>Polopriehľadný kruh</translation>
    </message>
    <message>
        <source>Cross</source>
        <translation>Krížik</translation>
    </message>
    <message>
        <source>Show all features</source>
        <translation type="unfinished">Ukázať všetky objekty</translation>
    </message>
    <message>
        <source>Show selected features</source>
        <translation type="unfinished">Ukázať vybrané objekty</translation>
    </message>
    <message>
        <source>Show features in current canvas</source>
        <translation type="unfinished">Ukázať objekty na aktuálnom plátne</translation>
    </message>
</context>
<context>
    <name>QgsOptionsBase</name>
    <message>
        <source>QGIS Options</source>
        <translation>QGIS Vlastnosti</translation>
    </message>
    <message>
        <source>Hide splash screen at startup</source>
        <translation>Pri štarte skryť úvodnú upútavku</translation>
    </message>
    <message>
        <source>&lt;b&gt;Note: &lt;/b&gt;Theme changes take effect the next time QGIS is started</source>
        <translation>&lt;b&gt;Poznámka: &lt;/b&gt;Zmena témy sa prejaví až po najbližšom spustení QGIS</translation>
    </message>
    <message>
        <source>&amp;Rendering</source>
        <translation>Vy&amp;kresľovanie</translation>
    </message>
    <message>
        <source>Map display will be updated (drawn) after this many features have been read from the data source</source>
        <translation>Mapový pohľad bude aktualizovaný (vykreslený) potom, čo takýto počet objektov bude načítaný zo zdroja údajov</translation>
    </message>
    <message>
        <source>Select Global Default ...</source>
        <translation>Nastaviť predvolené zobrazenie ...</translation>
    </message>
    <message>
        <source>Make lines appear less jagged at the expense of some drawing performance</source>
        <translation>Vyhladiť čiary na úrok nižšieho výkonu vykresľovania </translation>
    </message>
    <message>
        <source>Measure tool</source>
        <translation>Nástroj na meranie</translation>
    </message>
    <message>
        <source>Search radius</source>
        <translation>Polomer vyhľadávania</translation>
    </message>
    <message>
        <source>Pro&amp;jection</source>
        <translation type="obsolete">Mapové &amp;zobrazenie</translation>
    </message>
    <message>
        <source>When layer is loaded that has no projection information</source>
        <translation type="obsolete">Keď je nahrávaná vrstva bez informácie o zobrazení</translation>
    </message>
    <message>
        <source>By default new la&amp;yers added to the map should be displayed</source>
        <translation>Zobrazovať (vykresľovať) &amp;novopridané vrstvy do mapy</translation>
    </message>
    <message>
        <source>Fix problems with incorrectly filled polygons</source>
        <translation>Vyhnúť sa problému s nesprávne vypĺňanými polygónmi</translation>
    </message>
    <message>
        <source>Continuously redraw the map when dragging the legend/map divider</source>
        <translation>Pre posúvaní oddeľovača legendy a mapy priebežne prekresľovať mapu</translation>
    </message>
    <message>
        <source>&amp;Map tools</source>
        <translation>&amp;Mapové nástroje</translation>
    </message>
    <message>
        <source>%</source>
        <translation>%</translation>
    </message>
    <message>
        <source>Panning and zooming</source>
        <translation>Posun a zmena pohľadu</translation>
    </message>
    <message>
        <source>Zoom</source>
        <translation>Približovanie/oddaľovanie</translation>
    </message>
    <message>
        <source>Zoom and recenter</source>
        <translation>Približovanie/oddaľovanie a vycentrovanie</translation>
    </message>
    <message>
        <source>Nothing</source>
        <translation>Žiadna</translation>
    </message>
    <message>
        <source>&amp;General</source>
        <translation>&amp;Všeobecné</translation>
    </message>
    <message>
        <source>Locale</source>
        <translation type="unfinished">Regionálne nastavenie</translation>
    </message>
    <message>
        <source>Locale to use instead</source>
        <translation type="unfinished">namiesto toho použiť</translation>
    </message>
    <message>
        <source>Additional Info</source>
        <translation>Doplňujúce informácie</translation>
    </message>
    <message>
        <source>Detected active locale on your system:</source>
        <translation type="unfinished">Zistené aktívne regionálne nastavenie vásho systému:</translation>
    </message>
    <message>
        <source>Selecting this will unselect the &apos;make lines less&apos; jagged toggle</source>
        <translation type="unfinished">Výberom tohto zašrtávacieho políča sa zruší výber voľby &apos;Vyhladiť čiary&apos;</translation>
    </message>
    <message>
        <source>Digitizing</source>
        <translation>Digitalizácia</translation>
    </message>
    <message>
        <source>Rubberband</source>
        <translation type="unfinished">Vyberací obdĺžnik</translation>
    </message>
    <message>
        <source>Line width in pixels</source>
        <translation type="unfinished">Hrúbka čiary v pixeloch</translation>
    </message>
    <message>
        <source>Snapping</source>
        <translation type="unfinished">Zameriavanie</translation>
    </message>
    <message>
        <source>Zoom to mouse cursor</source>
        <translation type="unfinished">Približovať ku kurzoru myši</translation>
    </message>
    <message>
        <source>Project files</source>
        <translation>Súbory projektu</translation>
    </message>
    <message>
        <source>Prompt to save project changes when required</source>
        <translation type="unfinished">Vyzvať na uloženie súboru projektu, pokiaľ je to potrebné  </translation>
    </message>
    <message>
        <source>Warn when opening a project file saved with an older version of QGIS</source>
        <translation>Upozorniť pri otváraní súboru projektu uloženého v staršej verzii QGISu</translation>
    </message>
    <message>
        <source>Default Map Appearance (overridden by project properties)</source>
        <translation type="unfinished">Predvolený vzhľad mapy (pred týmto nastavením majú prednosť nastavenia vo Vlastnostiach projektu)</translation>
    </message>
    <message>
        <source>Selection color</source>
        <translation>Farba výberu</translation>
    </message>
    <message>
        <source>Background color</source>
        <translation>Farba pozadia</translation>
    </message>
    <message>
        <source>&amp;Application</source>
        <translation>&amp;Aplikácia</translation>
    </message>
    <message>
        <source>Icon theme</source>
        <translation>Téma ikon</translation>
    </message>
    <message>
        <source>Capitalise layer names in legend</source>
        <translation>Písať mená vrstiev v legende veľkými písmenami</translation>
    </message>
    <message>
        <source>Rendering behavior</source>
        <translation>Nastavenie prekresľovania</translation>
    </message>
    <message>
        <source>Number of features to draw before updating the display</source>
        <translation>Počet objektov ktoré sa vykreslia pred aktualizáciou zobrazenia</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Lucida Grande&apos;; font-size:13pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;span style=&quot; font-weight:600;&quot;&gt;Note:&lt;/span&gt; Use zero to prevent display updates until all features have been rendered&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Lucida Grande&apos;; font-size:13pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;span style=&quot; font-weight:600;&quot;&gt;Note:&lt;/span&gt; Zadanie nuly spôsobí, že zobrazenie nebude aktulizované kým nebudú vykreslené všetky objekty&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>Rendering quality</source>
        <translation>Kvalita vykresľovania</translation>
    </message>
    <message>
        <source>Zoom factor</source>
        <translation>Faktor zväčšenia</translation>
    </message>
    <message>
        <source>Mouse wheel action</source>
        <translation>Činnosť kolieska myši</translation>
    </message>
    <message>
        <source>Rubberband color</source>
        <translation type="unfinished">Farba označovacieho obdĺžnika</translation>
    </message>
    <message>
        <source>Ellipsoid for distance calculations</source>
        <translation>Eliposoid používaný pri výpočtoch vzdialenosti</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Lucida Grande&apos;; font-size:13pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;span style=&quot; font-weight:600;&quot;&gt;Note:&lt;/span&gt; Specify the search radius as a percentage of the map width&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Lucida Grande&apos;; font-size:13pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;span style=&quot; font-weight:600;&quot;&gt;Poznámka:&lt;/span&gt; Určuje polomer vyhľadávania v percentách šírky mapy&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>Search radius for identifying features and displaying map tips</source>
        <translation>Polomer vyhľadávania pri identifikácii objektov a zobrazovanie mapových tipov</translation>
    </message>
    <message>
        <source>Line width</source>
        <translation type="unfinished">Hrúbka čiary</translation>
    </message>
    <message>
        <source>Line colour</source>
        <translation type="unfinished">Farba línie</translation>
    </message>
    <message>
        <source>Default snap mode</source>
        <translation type="unfinished">Predvolený mód zameriavania</translation>
    </message>
    <message>
        <source>Default snapping tolerance in layer units</source>
        <translation type="unfinished">Predvolená tolerancia zameriavania v mapových jednotkách vrstvy</translation>
    </message>
    <message>
        <source>Search radius for vertex edits in layer units</source>
        <translation type="unfinished">Polomer zameriavania pre úpravu uzlov v mapových jedotkách</translation>
    </message>
    <message>
        <source>Vertex markers</source>
        <translation type="unfinished">Značky uzlov</translation>
    </message>
    <message>
        <source>Marker style</source>
        <translation type="unfinished">Štýl značky</translation>
    </message>
    <message>
        <source>Prompt for projection</source>
        <translation type="obsolete">Spýtať sa na mapové zobrazenie</translation>
    </message>
    <message>
        <source>Project wide default projection will be used</source>
        <translation type="obsolete">Použiť mapové zobrazenie projektu</translation>
    </message>
    <message>
        <source>Global default projection displa&amp;yed below will be used</source>
        <translation type="obsolete">Použiť nižšie nastavené &amp;všeobecné mapové zobrazenie</translation>
    </message>
    <message>
        <source>Override system locale</source>
        <translation type="unfinished">Nebrať do úvahy regionálne nastavenie systému</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Lucida Grande&apos;; font-size:13pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;span style=&quot; font-weight:600;&quot;&gt;Note:&lt;/span&gt; Enabling / changing overide on local requires an application restart&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Lucida Grande&apos;; font-size:13pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;&lt;span style=&quot; font-weight:600;&quot;&gt;Note:&lt;/span&gt; Zapnutie / zmena overide on local requires an application restart&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>Display classification attribute names in legend</source>
        <translation type="unfinished">Zobraziť meno atribútu použitého pre klasifikáciu v legende</translation>
    </message>
    <message>
        <source>&lt;b&gt;Note:&lt;/b&gt; Use zero to prevent display updates until all features have been rendered</source>
        <translation type="unfinished">&lt;b&gt;Poznámka:&lt;/b&gt; Zadanie nuly spôsobí, že zobrazenie nebude aktulizované kým nebudú vykreslené všetky objekty</translation>
    </message>
    <message>
        <source>&lt;b&gt;Note:&lt;/b&gt; Specify the search radius as a percentage of the map width</source>
        <translation type="unfinished">&lt;b&gt;Poznámka&lt;/b&gt; Určuje polomer vyhľadávania v percentách šírky mapy</translation>
    </message>
    <message>
        <source>&lt;b&gt;Note:&lt;/b&gt; Enabling / changing overide on local requires an application restart</source>
        <translation type="unfinished">&lt;b&gt;Poznámka:&lt;/b&gt; Zapnutie / zmena nastavenia vyžaduje rreštart aplikácie</translation>
    </message>
    <message>
        <source>Proxy</source>
        <translation>Proxy server</translation>
    </message>
    <message>
        <source>Use proxy for web access</source>
        <translation type="unfinished">Požitie Proxy servera pre webový prístup</translation>
    </message>
    <message>
        <source>Host</source>
        <translation type="unfinished">Hostiteľ</translation>
    </message>
    <message>
        <source>Port</source>
        <translation>Port</translation>
    </message>
    <message>
        <source>User</source>
        <translation>Používateľ</translation>
    </message>
    <message>
        <source>Leave this blank if no proxy username / password are required</source>
        <translation type="unfinished">Pokiaľ nie potrebné užívateľské meno / heslo k proxy, nechajte toto pole prázdne</translation>
    </message>
    <message>
        <source>Password</source>
        <translation>Heslo</translation>
    </message>
    <message>
        <source>Open attribute table in a dock window</source>
        <translation type="unfinished">Otvoriť tabuľku atribútov v dokovanom okne</translation>
    </message>
    <message>
        <source>CRS</source>
        <translation type="unfinished">CRS</translation>
    </message>
    <message>
        <source>When layer is loaded that has no coordinate reference system (CRS)</source>
        <translation type="unfinished">Pri nahrávaní vrstvy bez informácie o súradnicovom referenčnom systéme (CRS)</translation>
    </message>
    <message>
        <source>Prompt for CRS</source>
        <translation type="unfinished">Spýtať sa na CRS</translation>
    </message>
    <message>
        <source>Project wide default CRS will be used</source>
        <translation type="unfinished">Použiť súradnicový systém (CRS) projektu</translation>
    </message>
    <message>
        <source>Global default CRS displa&amp;yed below will be used</source>
        <translation type="unfinished">Použiť nižšie nastavený &amp;všeobecný súradnicový systém (CRS)</translation>
    </message>
    <message>
        <source>Attribute table behaviour</source>
        <translation type="unfinished">Správanie tabuľky atribútov</translation>
    </message>
</context>
<context>
    <name>QgsPasteTransformationsBase</name>
    <message>
        <source>Paste Transformations</source>
        <translation>Vložiť transformácie</translation>
    </message>
    <message>
        <source>&lt;b&gt;Note: This function is not useful yet!&lt;/b&gt;</source>
        <translation>&lt;b&gt;Poznámka: Táto funkcia zatiaľ nie je použiteľná!&lt;/b&gt;</translation>
    </message>
    <message>
        <source>Source</source>
        <translation>Zdroj</translation>
    </message>
    <message>
        <source>Destination</source>
        <translation>Cieľ</translation>
    </message>
    <message>
        <source>&amp;Help</source>
        <translation>&amp;Pomocník</translation>
    </message>
    <message>
        <source>F1</source>
        <translation>F1</translation>
    </message>
    <message>
        <source>Add New Transfer</source>
        <translation>Pridať nový prenos</translation>
    </message>
    <message>
        <source>&amp;OK</source>
        <translation>&amp;OK</translation>
    </message>
    <message>
        <source>&amp;Cancel</source>
        <translation>&amp;Zrušiť</translation>
    </message>
</context>
<context>
    <name>QgsPatternDialogBase</name>
    <message>
        <source>Select a fill pattern</source>
        <translation type="obsolete">Vybrať vzorku výplne</translation>
    </message>
    <message>
        <source>Cancel</source>
        <translation type="obsolete">Zrušiť</translation>
    </message>
    <message>
        <source>Ok</source>
        <translation type="obsolete">OK</translation>
    </message>
    <message>
        <source>No Fill</source>
        <translation type="obsolete">Bez výplne</translation>
    </message>
</context>
<context>
    <name>QgsPgGeoprocessing</name>
    <message>
        <source>Buffer features in layer %1</source>
        <translation>Vytvoriť okolie (buffer) objektov vo vrstve %1</translation>
    </message>
    <message>
        <source>Error connecting to the database</source>
        <translation>Nastala chyba pri pripájaní sa k databáze</translation>
    </message>
    <message>
        <source>&amp;Geoprocessing</source>
        <translation>&amp;Geoprocessing</translation>
    </message>
    <message>
        <source>Unable to add geometry column</source>
        <translation>Nemožno pridať stĺpec s geometriou</translation>
    </message>
    <message>
        <source>Unable to add geometry column to the output table </source>
        <translation> Do výstupnej tabuľky nie je možné pridať stĺpec s geometriou</translation>
    </message>
    <message>
        <source>Unable to create table</source>
        <translation>Nemožno vytvoriť tabuľku</translation>
    </message>
    <message>
        <source>Failed to create the output table </source>
        <translation> Pokus o vytvorenie výslednej tabuľky zlyhal</translation>
    </message>
    <message>
        <source>No GEOS support</source>
        <translation>Bez podpory GEOS-u</translation>
    </message>
    <message>
        <source>Buffer function requires GEOS support in PostGIS</source>
        <translation>Funkcia buffer (tvorba okolia) vyžaduje podporu GEOS-u v PostGISe </translation>
    </message>
    <message>
        <source>No Active Layer</source>
        <translation>Žiadna vrstva nie je aktívna</translation>
    </message>
    <message>
        <source>You must select a layer in the legend to buffer</source>
        <translation>Musíte vybrať vrstvu v legende pre ktorú budú vytvorené okolia (buffer)</translation>
    </message>
    <message>
        <source>A new layer is created in the database with the buffered features.</source>
        <translation>Nová vrstva je vytvorená v databáze s objektami s okolím (bufferom).</translation>
    </message>
    <message>
        <source>&amp;Buffer features</source>
        <translation>Vytvoriť &amp;okolie (buffer) objektov</translation>
    </message>
    <message>
        <source>Create a buffer for a PostgreSQL layer. </source>
        <translation type="unfinished">Vytvoriť okolie (buffer) pre vrstvu PostgreSQL. </translation>
    </message>
    <message>
        <source>Not a PostgreSQL/PostGIS Layer</source>
        <translation type="unfinished">Nie je vrstvou PostgreSQL/PostGIS</translation>
    </message>
    <message>
        <source> is not a PostgreSQL/PostGIS layer.
</source>
        <translation type="unfinished"> nie je vrstvou PostgreSQL/PostGIS.
</translation>
    </message>
    <message>
        <source>Geoprocessing functions are only available for PostgreSQL/PostGIS Layers</source>
        <translation type="unfinished">Funkcie geoprocessing-u sú dostupné len pre vrstvy PostgreSQL/PostGIS</translation>
    </message>
</context>
<context>
    <name>QgsPgQueryBuilder</name>
    <message>
        <source>Table &lt;b&gt;%1&lt;/b&gt; in database &lt;b&gt;%2&lt;/b&gt; on host &lt;b&gt;%3&lt;/b&gt;, user &lt;b&gt;%4&lt;/b&gt;</source>
        <translation>Tabuľka &lt;b&gt;%1&lt;/b&gt; v databáze &lt;b&gt;%2&lt;/b&gt; na hostiteľovi &lt;b&gt;%3&lt;/b&gt;, používateľ &lt;b&gt;%4&lt;/b&gt;</translation>
    </message>
    <message>
        <source>Query Result</source>
        <translation>Výsledok dopytu</translation>
    </message>
    <message>
        <source>The where clause returned </source>
        <translation>Klauzula WHERE vrátila</translation>
    </message>
    <message>
        <source> rows.</source>
        <translation>riadkov.</translation>
    </message>
    <message>
        <source>Query Failed</source>
        <translation>Dopyt zlyhal</translation>
    </message>
    <message>
        <source>An error occurred when executing the query:</source>
        <translation>Pri vykonávaní dopytu nastala chyba:</translation>
    </message>
    <message>
        <source>Connection Failed</source>
        <translation>Spojenie zlyhalo</translation>
    </message>
    <message>
        <source>Connection to the database failed:</source>
        <translation>Spojenie k databáze zlyhalo:</translation>
    </message>
    <message>
        <source>Database error</source>
        <translation>Chyba databázy</translation>
    </message>
    <message>
        <source>No Records</source>
        <translation>Žiadne záznamy</translation>
    </message>
    <message>
        <source>The query you specified results in zero records being returned. Valid PostgreSQL layers must have at least one feature.</source>
        <translation>Výsledkom vami určeného dopytu sú nulové (žiadne) záznamy. Platné vrstvy PostgreSQL však musia mať aspoň jeden objekt.
</translation>
    </message>
    <message>
        <source>&lt;p&gt;Failed to get sample of field values using SQL:&lt;/p&gt;&lt;p&gt;</source>
        <translation type="unfinished">&lt;p&gt;Nepodarilo sa získať vzorku údajov z jednotlivých polí s použitím SQL dopytu:&lt;/p&gt;&lt;p&gt;</translation>
    </message>
    <message>
        <source>No Query</source>
        <translation>Žiadny dopyt</translation>
    </message>
    <message>
        <source>You must create a query before you can test it</source>
        <translation type="unfinished">Je potrebné najprv vytvoriť dopyt predtým než ho budete testovať</translation>
    </message>
    <message>
        <source>Error in Query</source>
        <translation>Chyba v dopyte</translation>
    </message>
</context>
<context>
    <name>QgsPgQueryBuilderBase</name>
    <message>
        <source>PostgreSQL Query Builder</source>
        <translation>Nástroj na tvorbu PostreSQL dopytov</translation>
    </message>
    <message>
        <source>Clear</source>
        <translation>Vyčistiť</translation>
    </message>
    <message>
        <source>Test</source>
        <translation>Skúška</translation>
    </message>
    <message>
        <source>Ok</source>
        <translation>OK</translation>
    </message>
    <message>
        <source>Cancel</source>
        <translation>Zrušiť</translation>
    </message>
    <message>
        <source>SQL where clause</source>
        <translation>SQL klauzula WHERE</translation>
    </message>
    <message>
        <source>Operators</source>
        <translation>Operátory</translation>
    </message>
    <message>
        <source>=</source>
        <translation>=</translation>
    </message>
    <message>
        <source>IN</source>
        <translation>JE V</translation>
    </message>
    <message>
        <source>NOT IN</source>
        <translation>NIE JE V</translation>
    </message>
    <message>
        <source>&lt;</source>
        <translation>&lt;</translation>
    </message>
    <message>
        <source>&gt;</source>
        <translation>&gt;</translation>
    </message>
    <message>
        <source>%</source>
        <translation>%</translation>
    </message>
    <message>
        <source>&lt;=</source>
        <translation>&lt;=</translation>
    </message>
    <message>
        <source>&gt;=</source>
        <translation>&gt;=</translation>
    </message>
    <message>
        <source>!=</source>
        <translation>!=</translation>
    </message>
    <message>
        <source>LIKE</source>
        <translation>AKO</translation>
    </message>
    <message>
        <source>AND</source>
        <translation>A</translation>
    </message>
    <message>
        <source>ILIKE</source>
        <translation>NIE AKO</translation>
    </message>
    <message>
        <source>OR</source>
        <translation>ALEBO</translation>
    </message>
    <message>
        <source>NOT</source>
        <translation>NIE JE</translation>
    </message>
    <message>
        <source>Values</source>
        <translation>Hodnoty</translation>
    </message>
    <message>
        <source>All</source>
        <translation>Všetko</translation>
    </message>
    <message>
        <source>Sample</source>
        <translation>Vzorka</translation>
    </message>
    <message>
        <source>Fields</source>
        <translation>Polia</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:9pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Retrieve &lt;span style=&quot; font-weight:600;&quot;&gt;all&lt;/span&gt; the record in the vector file (&lt;span style=&quot; font-style:italic;&quot;&gt;if the table is big, the operation can consume some time&lt;/span&gt;)&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="unfinished">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:9pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Získať &lt;span style=&quot; font-weight:600;&quot;&gt;všetky&lt;/span&gt; záznamy z vektorového súboru (&lt;span style=&quot; font-style:italic;&quot;&gt;pokiaľ je tabuľka veľká, operácia môže trvať aj dlhší čas&lt;/span&gt;)&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:9pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Take a &lt;span style=&quot; font-weight:600;&quot;&gt;sample&lt;/span&gt; of records in the vector file&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="unfinished">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:9pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Získať &lt;span style=&quot; font-weight:600;&quot;&gt;vzorku&lt;/span&gt; záznamov z vektorového súboru&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:9pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;List of values for the current field.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="unfinished">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:9pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Zoznam hodnôt aktuálneho poľa.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:9pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;List of fields in this vector file&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="unfinished">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:9pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Zoznam polí v tomto vektorovom súbore&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>Datasource</source>
        <translation>Zdroj údajov</translation>
    </message>
</context>
<context>
    <name>QgsPluginInstaller</name>
    <message>
        <source>Couldn&apos;t parse output from the repository</source>
        <translation>Nemožno vyhodnotiť výstup z repozitára</translation>
    </message>
    <message>
        <source>Couldn&apos;t open the system plugin directory</source>
        <translation>Nemožno otvoriť adresár systémových zásuvných modulov</translation>
    </message>
    <message>
        <source>Couldn&apos;t open the local plugin directory</source>
        <translation>Nemožno otvoriť adresár lokálnych zásuvných modulov</translation>
    </message>
    <message>
        <source>Fetch Python Plugins...</source>
        <translation>Získať zásuvné moduly v Pythone...</translation>
    </message>
    <message>
        <source>Install more plugins from remote repositories</source>
        <translation>Nainštalovať viac zásuvných modulov zo vzdialených repozitárov</translation>
    </message>
    <message>
        <source>Looking for new plugins...</source>
        <translation>Hľadajú sa nové zásuvné moduly...</translation>
    </message>
    <message>
        <source>There is a new plugin available</source>
        <translation>Je dostupný nový zásuvný modul</translation>
    </message>
    <message>
        <source>There is a plugin update available</source>
        <translation>Je dostupná aktualizácia zásuvného modulu</translation>
    </message>
    <message>
        <source>QGIS Python Plugin Installer</source>
        <translation>QGIS Inštalátor zásuvných modulov v Pythone</translation>
    </message>
    <message>
        <source>Error reading repository:</source>
        <translation>Chyba pri čítaní repozitára:</translation>
    </message>
    <message>
        <source>Plugin directory doesn&apos;t exist:</source>
        <translation type="obsolete">Adresár so zásuvnými modulmi neexistuje:</translation>
    </message>
    <message>
        <source>Failed to remove the directory:</source>
        <translation>Nemožno odobrať adresár:</translation>
    </message>
    <message>
        <source>Check permissions or remove it manually</source>
        <translation>Skontrolujte práva alebo ho odoberte ručne</translation>
    </message>
    <message>
        <source>Nothing to remove! Plugin directory doesn&apos;t exist:</source>
        <translation>Niet čo odstrániť! Adresár so zásuvnými modulmi neexistuje:</translation>
    </message>
</context>
<context>
    <name>QgsPluginInstallerDialog</name>
    <message>
        <source>QGIS Python Plugin Installer</source>
        <translation>QGIS Inštalátor zásuvných modulov v Pythone</translation>
    </message>
    <message>
        <source>QGIS Plugin Installer</source>
        <translation type="obsolete">QGIS Inštalátor zásuvných modulov</translation>
    </message>
    <message>
        <source>Plugins</source>
        <translation type="obsolete">Zásuvné moduly</translation>
    </message>
    <message>
        <source>List of available and installed plugins</source>
        <translation type="obsolete">Zoznam dostupných a nainštalovaných zásuvných modulov</translation>
    </message>
    <message>
        <source>Filter:</source>
        <translation type="obsolete">Filter:</translation>
    </message>
    <message>
        <source>Display only plugins containing this word in their metadata</source>
        <translation type="obsolete">Zobraziť len zásuvné moduly obsahujúce toto slovo v svojich metaúdajoch</translation>
    </message>
    <message>
        <source>Display only plugins from given repository</source>
        <translation type="obsolete">Zobraziť len zásuvný modul z daného repozitára</translation>
    </message>
    <message>
        <source>all repositories</source>
        <translation>všetky repozitáre</translation>
    </message>
    <message>
        <source>Display only plugins with matching status</source>
        <translation type="obsolete">Zobraziť len zásuvné moduly so zodpovedajúcim stavom</translation>
    </message>
    <message>
        <source>Status</source>
        <translation type="obsolete">Stav</translation>
    </message>
    <message>
        <source>Name</source>
        <translation type="obsolete">Meno</translation>
    </message>
    <message>
        <source>Version</source>
        <translation type="obsolete">Verzia</translation>
    </message>
    <message>
        <source>Description</source>
        <translation type="obsolete">Popis</translation>
    </message>
    <message>
        <source>Author</source>
        <translation type="obsolete">Autor</translation>
    </message>
    <message>
        <source>Repository</source>
        <translation type="obsolete">Repozitár</translation>
    </message>
    <message>
        <source>Install, reinstall or upgrade the selected plugin</source>
        <translation type="obsolete">Inštalovať, reinštalovať alebo aktualizovať vybraný zásuvný modul</translation>
    </message>
    <message>
        <source>Install/upgrade plugin</source>
        <translation>Inštalovať/akualizovať zásuvný modul</translation>
    </message>
    <message>
        <source>Uninstall the selected plugin</source>
        <translation type="obsolete">Odinštalovať vybraný zásuvný modul</translation>
    </message>
    <message>
        <source>Uninstall plugin</source>
        <translation type="obsolete">Odinštalovať zásuvný modul</translation>
    </message>
    <message>
        <source>Repositories</source>
        <translation type="obsolete">Repozitáre</translation>
    </message>
    <message>
        <source>List of plugin repositories</source>
        <translation type="obsolete">Zoznam repozitárov zásuvných modulov</translation>
    </message>
    <message>
        <source>URL</source>
        <translation type="obsolete">URL</translation>
    </message>
    <message>
        <source>Allow the Installer to look for updates and news in enabled repositories on QGIS startup</source>
        <translation type="obsolete">Dovoliť Inštalátoru, aby pri štarte QGISu zisťoval aktualizácie a novinky v zapnutých repozitároch </translation>
    </message>
    <message>
        <source>Check for updates on startup</source>
        <translation type="obsolete">Pri štarte kontrolovať aktualizácie</translation>
    </message>
    <message>
        <source>Add third party plugin repositories to the list</source>
        <translation type="obsolete">Pridať do zoznamu repozitáre zásuvných modulov tretích strán</translation>
    </message>
    <message>
        <source>Add 3rd party repositories</source>
        <translation type="obsolete">Pridať repozitáre 3-ích strán</translation>
    </message>
    <message>
        <source>Add a new plugin repository</source>
        <translation type="obsolete">Pridá nový repozitár zásuvných modulov</translation>
    </message>
    <message>
        <source>Add...</source>
        <translation type="obsolete">Pridať...</translation>
    </message>
    <message>
        <source>Edit the selected repository</source>
        <translation type="obsolete">Upraví vybraný repozitár</translation>
    </message>
    <message>
        <source>Edit...</source>
        <translation type="obsolete">Upraviť...</translation>
    </message>
    <message>
        <source>Remove the selected repository</source>
        <translation type="obsolete">Odoberie vybraný repozitár</translation>
    </message>
    <message>
        <source>Delete</source>
        <translation type="obsolete">Vymazať</translation>
    </message>
    <message>
        <source>The plugins will be installed to ~/.qgis/python/plugins</source>
        <translation type="obsolete">Zásuvné moduly budú inštalované do ~/.qgis/python/plugins</translation>
    </message>
    <message>
        <source>Close the Installer window</source>
        <translation type="obsolete">Zatvorí okno inštalátora</translation>
    </message>
    <message>
        <source>Close</source>
        <translation type="obsolete">Zatvoriť</translation>
    </message>
    <message>
        <source>Error reading repository:</source>
        <translation>Chyba pri čítaní repozitára:</translation>
    </message>
    <message>
        <source>connected</source>
        <translation>pripojený</translation>
    </message>
    <message>
        <source>This repository is connected</source>
        <translation>Repozitár je pripojený</translation>
    </message>
    <message>
        <source>unavailable</source>
        <translation>nedostupný</translation>
    </message>
    <message>
        <source>This repository is enabled, but unavailable</source>
        <translation>Repozitár je zapnutý avšak je nedostupný</translation>
    </message>
    <message>
        <source>disabled</source>
        <translation>vypnutý</translation>
    </message>
    <message>
        <source>This repository is disabled</source>
        <translation>Repozitár je vypnutý</translation>
    </message>
    <message>
        <source>This repository is blocked due to incompatibility with your Quantum GIS version</source>
        <translation>Repozitár je zablokovaný kvôli nekompatibilite s vašou verziou QGISu</translation>
    </message>
    <message>
        <source>orphans</source>
        <translation>siroty</translation>
    </message>
    <message>
        <source>any status</source>
        <translation>akýkoľvek stav</translation>
    </message>
    <message>
        <source>not installed</source>
        <comment>plural</comment>
        <translation>nenainštalovaný</translation>
    </message>
    <message>
        <source>installed</source>
        <comment>plural</comment>
        <translation>nainštalovaný</translation>
    </message>
    <message>
        <source>upgradeable and news</source>
        <translation>aktualizácie a novinky</translation>
    </message>
    <message>
        <source>This plugin is not installed</source>
        <translation>Tento zásuvný modul nie je nainštalovaný</translation>
    </message>
    <message>
        <source>This plugin is installed</source>
        <translation>Tento zásuvný modul je nainštalovaný</translation>
    </message>
    <message>
        <source>This plugin is installed, but there is an updated version available</source>
        <translation type="unfinished">Tento zásuvný modul je nainštalovaný a je dostupná jeho aktualizácia</translation>
    </message>
    <message>
        <source>This plugin is installed, but I can&apos;t find it in any enabled repository</source>
        <translation type="unfinished">Tento zásuvný modul je nainštalovaný ale nemožno ho nájsť nikde v zapnutých repozitároch</translation>
    </message>
    <message>
        <source>This plugin is not installed and is seen for the first time</source>
        <translation type="unfinished">Tento zásuvný modul nie je nainštalovaný a objavil sa tu po prvý krát</translation>
    </message>
    <message>
        <source>This plugin is installed and is newer than its version available in a repository</source>
        <translation type="unfinished">Tento zásuvný modul je nainštalovaný a je novší ako jeho verzia dostupná v repozitároch</translation>
    </message>
    <message>
        <source>This plugin seems to be invalid or have unfulfilled dependencies
It has been installed, but can&apos;t be loaded</source>
        <translation type="obsolete">Tento zásuvný modul vyzerá byť chybný alebo nemá splnené závislosti
Bol nainštalovaný ale nemôže byť nahratý</translation>
    </message>
    <message>
        <source>not installed</source>
        <comment>singular</comment>
        <translation>nenainštalovaný</translation>
    </message>
    <message>
        <source>installed</source>
        <comment>singular</comment>
        <translation>nainštalovaný</translation>
    </message>
    <message>
        <source>upgradeable</source>
        <comment>singular</comment>
        <translation>aktualizovateľný</translation>
    </message>
    <message>
        <source>new!</source>
        <comment>singular</comment>
        <translation>nový!</translation>
    </message>
    <message>
        <source>invalid</source>
        <comment>singular</comment>
        <translation>chybný</translation>
    </message>
    <message>
        <source>installed version</source>
        <translation>verzia nainštalovaného</translation>
    </message>
    <message>
        <source>available version</source>
        <translation>dostupná verzia</translation>
    </message>
    <message>
        <source>That&apos;s the newest available version</source>
        <translation>To je najnovšia dostupná verzia</translation>
    </message>
    <message>
        <source>There is no version available for download</source>
        <translation>Na stiahnutie nie je dostupná žiadna verzia</translation>
    </message>
    <message>
        <source>only locally available</source>
        <translation>dostupné len lokálne</translation>
    </message>
    <message>
        <source>Install plugin</source>
        <translation>Nainštalovať zásuvný modul</translation>
    </message>
    <message>
        <source>Reinstall plugin</source>
        <translation>Preinštalovať zásuvný modul</translation>
    </message>
    <message>
        <source>Upgrade plugin</source>
        <translation>Aktualizovať zásuvný modul</translation>
    </message>
    <message>
        <source>Downgrade plugin</source>
        <translation>Deaktualizovať zásuvný modul</translation>
    </message>
    <message>
        <source>Plugin installation failed</source>
        <translation>Inštalácia zásuvného modulu zlyhala</translation>
    </message>
    <message>
        <source>Plugin installed successfully</source>
        <translation>Zásuvný modul bol úspešne nainštalovaný</translation>
    </message>
    <message>
        <source>Python plugin installed.
You have to enable it in the Plugin Manager.</source>
        <translation type="obsolete">Zásuvný modul bol nainštalovaný.
Môžno ho zapnúť v Správcovi zásuvných modulov.</translation>
    </message>
    <message>
        <source>Python plugin reinstalled.
You have to restart Quantum GIS to reload it.</source>
        <translation type="obsolete">Zásuvný modul bol preinštalovaný.
Je potrebné znovu sputiť Quantum GIS aby sa načítali.</translation>
    </message>
    <message>
        <source>Plugin uninstall failed</source>
        <translation>Odinštalovanie zásuvného modulu zlyhalo</translation>
    </message>
    <message>
        <source>Are you sure you want to uninstall the following plugin?</source>
        <translation type="unfinished">Ste si istý že chcete odinštalovať nasledujúci modul?</translation>
    </message>
    <message>
        <source>Warning: this plugin isn&apos;t available in any accessible repository!</source>
        <translation type="unfinished">Upozornenie: Tento zásuvný modul nie je dispozícii v žiadnom repozitári!</translation>
    </message>
    <message>
        <source>Plugin uninstalled successfully</source>
        <translation type="unfinished">Zásuvný modul bol úspešne odinštalovaný</translation>
    </message>
    <message>
        <source>You are going to add some plugin repositories neither authorized nor supported by the Quantum GIS team, however provided by folks associated with us. Plugin authors generally make efforts to make their works useful and safe, but we can&apos;t assume any responsibility for them. FEEL WARNED!</source>
        <translation type="obsolete">Chystáte sa pridať repozitáre zásuvných modulov ktoré nie sú autorizované ani podporvané tímom Quantum GIS ale ľuďmi s nami spojenými. Autori zásuvných modulov sa vo všeobecnosti snažia aby ich práce boli spoľahlivé a bezpečné ale nemôžeme niesť za nich žiadnu zodpovednosť. Takže boli ste varovaný!</translation>
    </message>
    <message>
        <source>Unable to add another repository with the same URL!</source>
        <translation type="unfinished">Nemožno pridať ďalší repozitár s rovankým URL!</translation>
    </message>
    <message>
        <source>Are you sure you want to remove the following repository?</source>
        <translation type="unfinished">Ste si istý, že chcete odobrať nasledujúci repozitár?</translation>
    </message>
    <message>
        <source>This plugin is incompatible with your Quantum GIS version and probably won&apos;t work.</source>
        <translation type="unfinished">Tento zásvuný modul je nekompatibilný s vašou verziou Quantum GISu a pravdepodobne nebude fungovať.</translation>
    </message>
    <message>
        <source>The required Python module is not installed.
For more information, please visit its homepage.</source>
        <translation type="obsolete">Požadovaný modul Pythonu nie je nainštalovaný.
Ďalšie informácie nájdete na jeho webovej stránke.</translation>
    </message>
    <message>
        <source>This plugin seems to be broken.
It has been installed but can&apos;t be loaded.
Here is the error message:</source>
        <translation type="unfinished">Tento zásuvný modul vyzerá byť poškodený.
Bol nainštalovaný ale nemožno ho nahrať.
Tu je chybové hlásenie:</translation>
    </message>
    <message>
        <source>Note that it&apos;s an uninstallable core plugin</source>
        <translation type="unfinished">Toto je neodinštalovateľný jadrový zásuvný modul</translation>
    </message>
    <message>
        <source>This plugin is broken</source>
        <translation type="unfinished">Zásuvný modul je poškodený</translation>
    </message>
    <message>
        <source>This plugin requires a newer version of Quantum GIS</source>
        <translation type="unfinished">Tento zásuvný modul vyžaduje novšiu verziu Quantum GISu</translation>
    </message>
    <message>
        <source>This plugin requires a missing module</source>
        <translation type="unfinished">Tento zásuvný modul vyžaduje modul ktorý chýba</translation>
    </message>
    <message>
        <source>Are you sure you want to downgrade the plugin to the latest available version? The installed one is newer!</source>
        <translation type="unfinished">Ste si istý že chcete prejsť na posledne dostupnú verziu zásuvného modulu? Nainštalovaná verzia je novšia!</translation>
    </message>
    <message>
        <source>Plugin has disappeared</source>
        <translation type="unfinished">Zásuvný modul zmizol</translation>
    </message>
    <message>
        <source>The plugin seems to have been installed but I don&apos;t know where. Probably the plugin package contained a wrong named directory.
Please search the list of installed plugins. I&apos;m nearly sure you&apos;ll find the plugin there, but I just can&apos;t determine which of them it is. It also means that I won&apos;t be able to determine if this plugin is installed and inform you about available updates. However the plugin may work. Please contact the plugin author and submit this issue.</source>
        <translation type="unfinished">Tento zásuvný modul bol nainštalovaný, ale nevedno kam. Balíček zásuvného modulu obsahoval pravdepodobne nesprávny názov adresára.
Posím pohľadajte zásuvný modul v zozname nainštalovaných modulov. Nie je isté, že ho tam nájdete, ale balíčkovací systém ho tam jendoducho nevie rozpoznať. To tiež znamená, že sa nedá rozpoznať či je daný zásuvný modul nainšatolvaný a ani informovať od jeho aktualizáciách. Avšak zásuvný modul môže fungvať. Prosím kontaktujte jeho autora a oznámte mu tento problém.</translation>
    </message>
    <message>
        <source>Plugin reinstalled successfully</source>
        <translation type="unfinished">Zásuvný modul bol úspešne preinštalovaný</translation>
    </message>
    <message>
        <source>The plugin is designed for a newer version of Quantum GIS. The minimum required version is:</source>
        <translation type="unfinished">Tento zásuvný modul bol navrhnutý pre novšiu verziu Quantum GISu. Minimálna požadovaná verzia:</translation>
    </message>
    <message>
        <source>The plugin depends on some components missing on your system. You need to install the following Python module in order to enable it:</source>
        <translation type="unfinished">Tento zásuvný modul je závislý na niektorých komponentoch, ktoré chýbajú na vašom systéme. ABy mohli fungovať je potrebné nainštalovať nasledujúce zásuvné moduly:</translation>
    </message>
    <message>
        <source>The plugin is broken. Python said:</source>
        <translation type="unfinished">Zásuvný modul je pokazený. Výstup Pythonu:</translation>
    </message>
    <message>
        <source>The required Python module is not installed.
For more information, please visit its homepage and Quantum GIS wiki.</source>
        <translation type="unfinished">Požadovaný modul Pythonu nie je nainštalovaný.
Viac informácií nájdete na jeho domovskej stránke alebo Wiki projektu Quantum GIS.</translation>
    </message>
    <message>
        <source>Python plugin installed.
Now you need to enable it in Plugin Manager.</source>
        <translation type="unfinished">Zásuvný modul Pythonu nainštalovaný.
teraz je potrebné ho zapnúť v Správcovi zásuvných modulov.</translation>
    </message>
    <message>
        <source>Python plugin reinstalled.
You need to restart Quantum GIS in order to reload it.</source>
        <translation>Zásuvný modul Pythonu preinštalovaný.
Aby mohol byť znova nahratý, je potrebné najprv reštartovať Quantum GIS.</translation>
    </message>
    <message>
        <source>Python plugin uninstalled. Note that tou may need to restart Quantum GIS in order to remove it completely.</source>
        <translation type="obsolete">Zásuvný modul Pythonu odinštalovaný. Na to, aby bol úplne odinštalovaný je potrebné reštartovať Quantum GIS.</translation>
    </message>
    <message>
        <source>Python plugin uninstalled. Note that you may need to restart Quantum GIS in order to remove it completely.</source>
        <translation>Zásuvný modul v Pythone bol odinštalovaný. Na úplné dokončenie odinštalovania bude potrebné reštartovať Quantum GIS.</translation>
    </message>
    <message>
        <source>You are about to add several plugin repositories that are neither authorized nor supported by the Quantum GIS team. Plugin authors generally make efforts to ensure that their work is useful and safe, however, we can assume no responsibility for them.</source>
        <translation type="unfinished"></translation>
    </message>
</context>
<context>
    <name>QgsPluginInstallerDialogBase</name>
    <message>
        <source>QGIS Python Plugin Installer</source>
        <translation type="unfinished">QGIS Inštalátor zásuvných modulov Pythonu</translation>
    </message>
    <message>
        <source>Plugins</source>
        <translation>Zásuvné moduly</translation>
    </message>
    <message>
        <source>List of available and installed plugins</source>
        <translation>Zoznam dostupných a nainštalovaných zásuvných modulov</translation>
    </message>
    <message>
        <source>Filter:</source>
        <translation>Filter:</translation>
    </message>
    <message>
        <source>Display only plugins containing this word in their metadata</source>
        <translation type="unfinished">Zobraziť len zásuvné moduly obsahujúce toto slovo v svojich metaúdajoch</translation>
    </message>
    <message>
        <source>Display only plugins from given repository</source>
        <translation type="unfinished">Zobraziť len zásuvný modul z daného repozitára</translation>
    </message>
    <message>
        <source>all repositories</source>
        <translation>všetky repozitáre</translation>
    </message>
    <message>
        <source>Display only plugins with matching status</source>
        <translation type="unfinished">Zobraziť len zásuvné moduly so zodpovedajúcim stavom</translation>
    </message>
    <message>
        <source>Status</source>
        <translation>Stav</translation>
    </message>
    <message>
        <source>Name</source>
        <translation type="unfinished">Meno</translation>
    </message>
    <message>
        <source>Version</source>
        <translation>Verzia</translation>
    </message>
    <message>
        <source>Description</source>
        <translation>Popis</translation>
    </message>
    <message>
        <source>Author</source>
        <translation>Autor</translation>
    </message>
    <message>
        <source>Repository</source>
        <translation>Repozitár</translation>
    </message>
    <message>
        <source>Install, reinstall or upgrade the selected plugin</source>
        <translation>Nainštaluje, preinštaluje, alebo aktualizuje vybraný zásuvný modul</translation>
    </message>
    <message>
        <source>Install/upgrade plugin</source>
        <translation>Inštalovať/akualizovať zásuvný modul</translation>
    </message>
    <message>
        <source>Uninstall the selected plugin</source>
        <translation>Odinštaluje vybraný zásuvný modul</translation>
    </message>
    <message>
        <source>Uninstall plugin</source>
        <translation>Odinštalovať zásuvný modul</translation>
    </message>
    <message>
        <source>Repositories</source>
        <translation>Repozitáre</translation>
    </message>
    <message>
        <source>List of plugin repositories</source>
        <translation type="unfinished">Zoznam repozitárov zásuvných modulov</translation>
    </message>
    <message>
        <source>URL</source>
        <translation>URL</translation>
    </message>
    <message>
        <source>Allow the Installer to look for updates and news in enabled repositories on QGIS startup</source>
        <translation type="unfinished">Dovoliť Inštalátoru, aby pri štarte QGISu zisťoval aktualizácie a novinky v zapnutých repozitároch </translation>
    </message>
    <message>
        <source>Check for updates on startup</source>
        <translation>Pri štarte kontrolovať aktualizácie</translation>
    </message>
    <message>
        <source>Add third party plugin repositories to the list</source>
        <translation>Pridá do zoznamu repozitáre zásuvných modulov tretích strán</translation>
    </message>
    <message>
        <source>Add 3rd party repositories</source>
        <translation>Pridať repozitáre 3-ích strán</translation>
    </message>
    <message>
        <source>Add a new plugin repository</source>
        <translation>Pridá nový repozitár zásuvných modulov</translation>
    </message>
    <message>
        <source>Add...</source>
        <translation>Pridať...</translation>
    </message>
    <message>
        <source>Edit the selected repository</source>
        <translation>Upraví vybraný repozitár</translation>
    </message>
    <message>
        <source>Edit...</source>
        <translation>Upraviť...</translation>
    </message>
    <message>
        <source>Remove the selected repository</source>
        <translation>Odoberie vybraný repozitár</translation>
    </message>
    <message>
        <source>Delete</source>
        <translation>Vymazať</translation>
    </message>
    <message>
        <source>The plugins will be installed to ~/.qgis/python/plugins</source>
        <translation>Zásuvné moduly budú nainštalované do ~/.qgis/python/plugins</translation>
    </message>
    <message>
        <source>Close the Installer window</source>
        <translation>Zatvorí okno inštalátora</translation>
    </message>
    <message>
        <source>Close</source>
        <translation>Zatvoriť</translation>
    </message>
</context>
<context>
    <name>QgsPluginInstallerFetchingDialog</name>
    <message>
        <source>Fetching repositories</source>
        <translation type="obsolete">Sťahujú sa repozitáre</translation>
    </message>
    <message>
        <source>Overall progress:</source>
        <translation type="obsolete">Celkový priebeh:</translation>
    </message>
    <message>
        <source>Abort fetching</source>
        <translation type="obsolete">Zrušiť sťahovanie</translation>
    </message>
    <message>
        <source>Repository</source>
        <translation type="obsolete">Repozitár</translation>
    </message>
    <message>
        <source>State</source>
        <translation type="obsolete">Stav</translation>
    </message>
    <message>
        <source>Success</source>
        <translation type="unfinished">Úspech</translation>
    </message>
    <message>
        <source>Resolving host name...</source>
        <translation type="unfinished">Zisťuje sa meno servera...</translation>
    </message>
    <message>
        <source>Connecting...</source>
        <translation>Pripája sa...</translation>
    </message>
    <message>
        <source>Host connected. Sending request...</source>
        <translation type="unfinished">Spojený so serverom. Odosiela sa požiadavka...</translation>
    </message>
    <message>
        <source>Downloading data...</source>
        <translation type="unfinished">Sťahujú sa údaje...</translation>
    </message>
    <message>
        <source>Idle</source>
        <translation type="unfinished">Idle</translation>
    </message>
    <message>
        <source>Closing connection...</source>
        <translation type="unfinished">Zatvára sa spojenie...</translation>
    </message>
    <message>
        <source>Error</source>
        <translation>Chyba</translation>
    </message>
</context>
<context>
    <name>QgsPluginInstallerFetchingDialogBase</name>
    <message>
        <source>Fetching repositories</source>
        <translation type="unfinished">Sťahujú sa repozitáre</translation>
    </message>
    <message>
        <source>Overall progress:</source>
        <translation type="unfinished">Celkový priebeh:</translation>
    </message>
    <message>
        <source>Abort fetching</source>
        <translation type="unfinished">Zrušiť sťahovanie</translation>
    </message>
    <message>
        <source>Repository</source>
        <translation>Repozitár</translation>
    </message>
    <message>
        <source>State</source>
        <translation type="unfinished">Stav</translation>
    </message>
</context>
<context>
    <name>QgsPluginInstallerInstallingDialog</name>
    <message>
        <source>QGIS Python Plugin Installer</source>
        <translation type="obsolete">QGIS Inštalátor zásuvných modulov Pythonu</translation>
    </message>
    <message>
        <source>Installing plugin:</source>
        <translation type="obsolete">Inštalovať zásuvný modul:</translation>
    </message>
    <message>
        <source>Connecting...</source>
        <translation>Pripája sa...</translation>
    </message>
    <message>
        <source>Installing...</source>
        <translation>Inštaluje sa...</translation>
    </message>
    <message>
        <source>Resolving host name...</source>
        <translation>Zisťuje sa meno servera...</translation>
    </message>
    <message>
        <source>Host connected. Sending request...</source>
        <translation>Server pripojený. Odosiela sa požiadavka...</translation>
    </message>
    <message>
        <source>Downloading data...</source>
        <translation>Sťahujú sa údaje...</translation>
    </message>
    <message>
        <source>Idle</source>
        <translation>Nečinný</translation>
    </message>
    <message>
        <source>Closing connection...</source>
        <translation>Zatvára sa spojenie...</translation>
    </message>
    <message>
        <source>Error</source>
        <translation>Chyba</translation>
    </message>
    <message>
        <source>Failed to unzip file to the following directory:</source>
        <translation type="obsolete">Rozpakovaie súboru do nasledujúceho adresára zlyhalo:</translation>
    </message>
    <message>
        <source>Check permissions</source>
        <translation type="obsolete">Skontrolujte práva</translation>
    </message>
    <message>
        <source>Aborted by user</source>
        <translation>Zrušené používateľom</translation>
    </message>
    <message>
        <source>Failed to unzip the plugin package. Probably it&apos;s broken or missing from the repository. You may also want to make sure that you have write permission to the plugin directory:</source>
        <translation>Balík so zásuvným modulom sa nepodarilo rozbaliť. Pravdepodobne je pokazený alebo chýba vo vašom repozitári. Tiež je dobré uistiť sa, že v adresári zásuvných modulov je nastavené právo na zápis:</translation>
    </message>
</context>
<context>
    <name>QgsPluginInstallerInstallingDialogBase</name>
    <message>
        <source>QGIS Python Plugin Installer</source>
        <translation>QGIS Inštalátor zásuvných modulov Pythonu</translation>
    </message>
    <message>
        <source>Installing plugin:</source>
        <translation>Inštalácia zásuvného modulu:</translation>
    </message>
    <message>
        <source>Connecting...</source>
        <translation>Pripája sa...</translation>
    </message>
</context>
<context>
    <name>QgsPluginInstallerPluginErrorDialog</name>
    <message>
        <source>Error loading plugin</source>
        <translation type="obsolete">Chyba pri nahrávaní zásuvného modulu</translation>
    </message>
    <message>
        <source>The plugin seems to be invalid or have unfulfilled dependencies. It has been installed, but can&apos;t be loaded. If you really need this plugin, you can contact its author or &lt;a href=&quot;http://lists.osgeo.org/mailman/listinfo/qgis-user&quot;&gt;QGIS users group&lt;/a&gt; and try to solve the problem. If not, you can just uninstall it. Here is the error message below:</source>
        <translation type="obsolete">Tento zásuvný modul vyzerá byť chybný alebo má nesplnené závislosti. Bol nainštalovaný ale nemožno ho nahrať. Ak naozaj potrebujete tento zásuvný modul, môžete kontaktovať autora alebo &lt;a href=&quot;http://lists.osgeo.org/mailman/listinfo/qgis-user&quot;&gt;skupinu užívateľov QGISu&lt;/a&gt; a pokúsiť sa vyriešiť tento problém. V opačnom prípade ho stačí odinštalovať. Tu je chybové hlásenie:</translation>
    </message>
    <message>
        <source>Do you want to uninstall this plugin now? If you&apos;re unsure, probably you would like to do this.</source>
        <translation type="obsolete">Želáte si odinštalovať tento zásuvný modul? Ak si nie ste istý, pravdepodobne ste to chceli.</translation>
    </message>
    <message>
        <source>No error message received. Try to restart Quantum GIS and ensure the plugin isn&apos;t installed under a different name. If it is, contact the plugin author and submit this issue, please.</source>
        <translation type="obsolete">Nebola obdržaná žiadna chybová správa. Pokúste sa znovuspustiť Quantum GIS a uistiť sa že tento zásuvný modul nie je nainštalovaný pod iným menom. Pokiaľ je, kontaktujte autora zásuvného modulu a oznámte mu túto záležitosť.</translation>
    </message>
    <message>
        <source>no error message received</source>
        <translation type="unfinished">nebola prijatá žiadna chybová správa</translation>
    </message>
</context>
<context>
    <name>QgsPluginInstallerPluginErrorDialogBase</name>
    <message>
        <source>Error loading plugin</source>
        <translation>Chyba pri nahrávaní zásuvného modulu</translation>
    </message>
    <message>
        <source>The plugin seems to be invalid or have unfulfilled dependencies. It has been installed, but can&apos;t be loaded. If you really need this plugin, you can contact its author or &lt;a href=&quot;http://lists.osgeo.org/mailman/listinfo/qgis-user&quot;&gt;QGIS users group&lt;/a&gt; and try to solve the problem. If not, you can just uninstall it. Here is the error message below:</source>
        <translation type="unfinished">Tento zásuvný modul vyzerá byť chybný alebo má nesplnené závislosti. Bol nainštalovaný ale nemožno ho nahrať. Ak naozaj potrebujete tento zásuvný modul, môžete kontaktovať autora alebo &lt;a href=&quot;http://lists.osgeo.org/mailman/listinfo/qgis-user&quot;&gt;skupinu užívateľov QGISu&lt;/a&gt; a pokúsiť sa vyriešiť tento problém. V opačnom prípade ho stačí odinštalovať. Tu je chybové hlásenie:</translation>
    </message>
    <message>
        <source>Do you want to uninstall this plugin now? If you&apos;re unsure, probably you would like to do this.</source>
        <translation type="unfinished">Želáte si odinštalovať tento zásuvný modul? Ak si nie ste istý, pravdepodobne ste to chceli.</translation>
    </message>
</context>
<context>
    <name>QgsPluginInstallerRepositoryDetailsDialog</name>
    <message>
        <source>Repository details</source>
        <translation type="obsolete">Podrobnosti repozitára</translation>
    </message>
    <message>
        <source>Name:</source>
        <translation type="obsolete">Meno:</translation>
    </message>
    <message>
        <source>Enter a name for the repository</source>
        <translation type="obsolete">Zadajte meno repozitára</translation>
    </message>
    <message>
        <source>URL:</source>
        <translation type="obsolete">URL:</translation>
    </message>
    <message>
        <source>Enter the repository URL, beginning with &quot;http://&quot;</source>
        <translation type="obsolete">Zadajte URL repozitára začínajúce s &quot;http://&quot;</translation>
    </message>
    <message>
        <source>Enable or disable the repository (disabled repositories will be omitted)</source>
        <translation type="obsolete">Povoľte alebo zrušte povolenie pre repozitár (zrušené repozitáre budú poreskočené)</translation>
    </message>
    <message>
        <source>Enabled</source>
        <translation type="obsolete">Povolený</translation>
    </message>
    <message>
        <source>[place for a warning message]</source>
        <translation type="obsolete">[miesto pre upozornenia]</translation>
    </message>
</context>
<context>
    <name>QgsPluginInstallerRepositoryDetailsDialogBase</name>
    <message>
        <source>Repository details</source>
        <translation type="unfinished">Podrobnosti o repozitári</translation>
    </message>
    <message>
        <source>Name:</source>
        <translation type="unfinished">Meno:</translation>
    </message>
    <message>
        <source>Enter a name for the repository</source>
        <translation type="unfinished">Zadajte meno repozitára</translation>
    </message>
    <message>
        <source>URL:</source>
        <translation>URL:</translation>
    </message>
    <message>
        <source>Enter the repository URL, beginning with &quot;http://&quot;</source>
        <translation type="unfinished">Zadajte URL repozitára začínajúce s &quot;http://&quot;</translation>
    </message>
    <message>
        <source>Enable or disable the repository (disabled repositories will be omitted)</source>
        <translation type="unfinished">Povoľte alebo zrušte povolenie pre repozitár (zrušené repozitáre budú poreskočené)</translation>
    </message>
    <message>
        <source>Enabled</source>
        <translation type="unfinished">Povolený</translation>
    </message>
</context>
<context>
    <name>QgsPluginManager</name>
    <message>
        <source>No Plugins</source>
        <translation>Žiadne zásuvné moduly</translation>
    </message>
    <message>
        <source>No QGIS plugins found in </source>
        <translation>Nenašli sa žiadne zásuvné moduly adresári </translation>
    </message>
    <message>
        <source>&amp;Select All</source>
        <translation>&amp;Vybrať všetky</translation>
    </message>
    <message>
        <source>&amp;Clear All</source>
        <translation>&amp;Nevybrať žiadny</translation>
    </message>
    <message>
        <source>[ incompatible ]</source>
        <translation>[ nekompatibilný ]</translation>
    </message>
</context>
<context>
    <name>QgsPluginManagerBase</name>
    <message>
        <source>QGIS Plugin Manager</source>
        <translation>Správca zásuvných modulov QGISu</translation>
    </message>
    <message>
        <source>To enable / disable a plugin, click its checkbox or description</source>
        <translation>Kliknutím na zaškrtávacie políčko alebo popis možno zapnúť alebo vypnúť zásuvný modul</translation>
    </message>
    <message>
        <source>&amp;Filter</source>
        <translation>&amp;Filter</translation>
    </message>
    <message>
        <source>Plugin Directory:</source>
        <translation>Adresár zásuvných modulov:</translation>
    </message>
    <message>
        <source>Directory</source>
        <translation>Adresár</translation>
    </message>
</context>
<context>
    <name>QgsPointDialog</name>
    <message>
        <source>Zoom In</source>
        <translation>Priblížiť</translation>
    </message>
    <message>
        <source>z</source>
        <translation>z</translation>
    </message>
    <message>
        <source>Zoom Out</source>
        <translation>Oddialiť</translation>
    </message>
    <message>
        <source>Z</source>
        <translation>Z</translation>
    </message>
    <message>
        <source>Zoom To Layer</source>
        <translation>Pohľad na veľkosť vrstvy</translation>
    </message>
    <message>
        <source>Zoom to Layer</source>
        <translation>Zmeniť pohľad na veľkosť vrstvy</translation>
    </message>
    <message>
        <source>Pan Map</source>
        <translation>Posun mapy</translation>
    </message>
    <message>
        <source>Pan the map</source>
        <translation>Posunúť mapu</translation>
    </message>
    <message>
        <source>Add Point</source>
        <translation>Pridať bod</translation>
    </message>
    <message>
        <source>.</source>
        <translation>.</translation>
    </message>
    <message>
        <source>Capture Points</source>
        <translation>Získať body</translation>
    </message>
    <message>
        <source>Delete Point</source>
        <translation>Zmazať bod</translation>
    </message>
    <message>
        <source>Delete Selected</source>
        <translation>Zmazať vybrané</translation>
    </message>
    <message>
        <source>Linear</source>
        <translation>Lineárna</translation>
    </message>
    <message>
        <source>Helmert</source>
        <translation>Helmertova</translation>
    </message>
    <message>
        <source>Choose a name for the world file</source>
        <translation type="unfinished">Vyberte meno world súboru</translation>
    </message>
    <message>
        <source>-modified</source>
        <comment>

Georeferencer:QgsPointDialog.cpp - used to modify a user given filename</comment>
        <translation type="obsolete">-upraveny</translation>
    </message>
    <message>
        <source>Warning</source>
        <translation>Upozornenie</translation>
    </message>
    <message>
        <source>Affine</source>
        <translation>Afinná</translation>
    </message>
    <message>
        <source>Not implemented!</source>
        <translation>Neimplementované!</translation>
    </message>
    <message>
        <source>&lt;p&gt;An affine transform requires changing the original raster file. This is not yet supported.&lt;/p&gt;</source>
        <translation>&lt;p&gt;Afinná transformácia vyžaduje zmeniť pôvodný rastrový súbor. Táto funkcia zatiaľ nie je implementovaná.&lt;/p&gt;</translation>
    </message>
    <message>
        <source>&lt;p&gt;The </source>
        <translation>&lt;p&gt;</translation>
    </message>
    <message>
        <source> transform is not yet supported.&lt;/p&gt;</source>
        <translation> transformácia nie je zatiaľ podporovaná.&lt;/p&gt;</translation>
    </message>
    <message>
        <source>Error</source>
        <translation>Chyba</translation>
    </message>
    <message>
        <source>Could not write to </source>
        <translation>Nemožno zapisovať do </translation>
    </message>
    <message>
        <source>Currently all modified files will be written in TIFF format.</source>
        <translation type="unfinished">V súčastnosti všetky zmenené súbory budú uložené vo formáte TIFF.</translation>
    </message>
    <message>
        <source>&lt;p&gt;A Helmert transform requires modifications in the raster layer.&lt;/p&gt;&lt;p&gt;The modified raster will be saved in a new file and a world file will be generated for this new file instead.&lt;/p&gt;&lt;p&gt;Are you sure that this is what you want?&lt;/p&gt;</source>
        <translation type="unfinished">&lt;p&gt;Helmertova transformácia vyžaduje úpravy v tejto rastrovej vrstve.&lt;/p&gt;&lt;p&gt;Pozmenený raster bude uložený v novom súbore a world súbor bude vygenerovaný pre tento nový súbor.&lt;/p&gt;&lt;p&gt;Ste si istý, že chcete vykonať túto operáciu?&lt;/p&gt;</translation>
    </message>
    <message>
        <source>-modified</source>
        <comment>Georeferencer:QgsPointDialog.cpp - used to modify a user given filename</comment>
        <translation type="obsolete">-upraveny</translation>
    </message>
    <message>
        <source>-modified</source>
        <comment>Georeferencer:QgsPointDialog.cpp - used to modify a user given file name</comment>
        <translation type="unfinished">-upravený</translation>
    </message>
</context>
<context>
    <name>QgsPointDialogBase</name>
    <message>
        <source>Transform type:</source>
        <translation>Typ transformácie:</translation>
    </message>
    <message>
        <source>Zoom in</source>
        <translation>Priblížiť</translation>
    </message>
    <message>
        <source>Zoom out</source>
        <translation>Oddialiť</translation>
    </message>
    <message>
        <source>Zoom to the raster extents</source>
        <translation>Zoom na rozsah rastra</translation>
    </message>
    <message>
        <source>Pan</source>
        <translation>Posun</translation>
    </message>
    <message>
        <source>Add points</source>
        <translation>Pridať body</translation>
    </message>
    <message>
        <source>Delete points</source>
        <translation>Vymazať body</translation>
    </message>
    <message>
        <source>World file:</source>
        <translation>World súbor:</translation>
    </message>
    <message>
        <source>...</source>
        <translation>...</translation>
    </message>
    <message>
        <source>Modified raster:</source>
        <translation>Pozmenený raster:</translation>
    </message>
    <message>
        <source>Reference points</source>
        <translation type="unfinished">Lícovacie body</translation>
    </message>
    <message>
        <source>Create</source>
        <translation type="unfinished">Vytvoriť</translation>
    </message>
    <message>
        <source>Create and load layer</source>
        <translation type="unfinished">Vytvoriť a nahrať vrstvu</translation>
    </message>
</context>
<context>
    <name>QgsPointStyleWidgetBase</name>
    <message>
        <source>Form3</source>
        <translation type="obsolete">Štýl bodu</translation>
    </message>
    <message>
        <source>Symbol Style</source>
        <translation type="obsolete">Štýl symbolu</translation>
    </message>
    <message>
        <source>Scale</source>
        <translation type="obsolete">Mierka</translation>
    </message>
</context>
<context>
    <name>QgsPostgresProvider</name>
    <message>
        <source>Unable to access relation</source>
        <translation>Nie je možné pristúpiť k relácii</translation>
    </message>
    <message>
        <source>Unable to access the </source>
        <translation>Nemožno pristúpiť k relácii </translation>
    </message>
    <message>
        <source> relation.
The error message from the database was:
</source>
        <translation> .
Chybové hlásenie z databázy bolo:
</translation>
    </message>
    <message>
        <source>No GEOS Support!</source>
        <translation type="obsolete">Bez podpory GEOS!</translation>
    </message>
    <message>
        <source>Your PostGIS installation has no GEOS support.
Feature selection and identification will not work properly.
Please install PostGIS with GEOS support (http://geos.refractions.net)</source>
        <translation type="obsolete">Vaša inštalácia PostGIS nemá podporu GEOSu.
Výber objektov a identifikácia nebudú pracovať správne.
Prosím nainštalujte PostGIS s podporou GEOSu (http://geos.refractions.net)</translation>
    </message>
    <message>
        <source>No suitable key column in table</source>
        <translation>V tabuľke nie je žiadny vhodný kľúčový stĺpec</translation>
    </message>
    <message>
        <source>The table has no column suitable for use as a key.

Qgis requires that the table either has a column of type
int4 with a unique constraint on it (which includes the
primary key) or has a PostgreSQL oid column.
</source>
        <translation>Táto tabuľka nemá žiadny stĺpec vhodný ako kľúč.

Qgis vyžaduje, aby tabuľka mala buď stĺpec typu
int4 s obmedzením na jedinečnosť (vrátane 
primárneho kľúča) alebo mala stĺpec PostreSQL oid.
</translation>
    </message>
    <message>
        <source>The unique index on column</source>
        <translation>Tento jedinečný index v stĺpci</translation>
    </message>
    <message>
        <source>is unsuitable because Qgis does not currently support non-int4 type columns as a key into the table.
</source>
        <translation>je nevhodný, pretože Qgis v súčasnosti nepodporuje ako kľúč tabuľky stĺpce iného typu ako int4.
 </translation>
    </message>
    <message>
        <source>and </source>
        <translation> a</translation>
    </message>
    <message>
        <source>The unique index based on columns </source>
        <translation> Tento jedinečný index založený na stĺpcoch</translation>
    </message>
    <message>
        <source> is unsuitable because Qgis does not currently support multiple columns as a key into the table.
</source>
        <translation> je nevhodný, pretože Qgis v súčasnosti nepodporuje ako kľúč tabuľky viacero stĺpcov.
 </translation>
    </message>
    <message>
        <source>Unable to find a key column</source>
        <translation>Nie je možné nájsť kľúčový stĺpec</translation>
    </message>
    <message>
        <source> derives from </source>
        <translation>  odvodený od</translation>
    </message>
    <message>
        <source>and is suitable.</source>
        <translation>a je vhodný.</translation>
    </message>
    <message>
        <source>and is not suitable </source>
        <translation> a nie je vhodný</translation>
    </message>
    <message>
        <source>type is </source>
        <translation>  typ je</translation>
    </message>
    <message>
        <source> and has a suitable constraint)</source>
        <translation> a má vhodné obmedzenie)</translation>
    </message>
    <message>
        <source> and does not have a suitable constraint)</source>
        <translation> a nemá vhodné obmedzenie)</translation>
    </message>
    <message>
        <source>The view you selected has the following columns, none of which satisfy the above conditions:</source>
        <translation>Pohľad, ktorý bol vybratý má nasledujúce stĺpce, žiadny z nich nevyhovuje vyššie uvedeným podmienkam:</translation>
    </message>
    <message>
        <source>Qgis requires that the view has a column that can be used as a unique key. Such a column should be derived from a table column of type int4 and be a primary key, have a unique constraint on it, or be a PostgreSQL oid column. To improve performance the column should also be indexed.
</source>
        <translation>Qgis vyžaduje, aby pohľad mal stĺpec, ktorý možno použiť ako jedinečný kľúč. Taký stĺpec by mal byť odvodený z tabuľky zo stĺpca typu int4 a byť primárnym kľúčom, mať obmedzenie jedinečnosti, alebo byť stĺpcom PostgreSQL oid. Pre zlepšenie výkonu by stĺpec mal byť indexovaný.
</translation>
    </message>
    <message>
        <source>The view </source>
        <translation>  Tento pohľad</translation>
    </message>
    <message>
        <source>has no column suitable for use as a unique key.
</source>
        <translation>nemá žiadny stĺpec použiteľný ako jedinečný kľúč.
</translation>
    </message>
    <message>
        <source>No suitable key column in view</source>
        <translation>Žiadny vhodný kľúčový stĺpec v aktuálnom pohľade</translation>
    </message>
    <message>
        <source>Unknown geometry type</source>
        <translation>Neznámy typ geometrie</translation>
    </message>
    <message>
        <source>Column </source>
        <translation> Stĺpec</translation>
    </message>
    <message>
        <source> in </source>
        <translation>  v</translation>
    </message>
    <message>
        <source> has a geometry type of </source>
        <translation>  obsahuje typ geometrie</translation>
    </message>
    <message>
        <source>, which Qgis does not currently support.</source>
        <translation>, s ktorým Qgis zatiaľ nevie pracovať.</translation>
    </message>
    <message>
        <source>. The database communication log was:
</source>
        <translation>. Správa z komunikácie s databázou:
</translation>
    </message>
    <message>
        <source>Unable to get feature type and srid</source>
        <translation>Nie je možné získať typ objektu a srid</translation>
    </message>
    <message>
        <source>Note: </source>
        <translation> Poznámka:</translation>
    </message>
    <message>
        <source>initially appeared suitable but does not contain unique data, so is not suitable.
</source>
        <translation>na začiatku sa javil ako vhodný, ale neobsahuje jedinečné údaje a preto vhodný nie je.
</translation>
    </message>
    <message>
        <source>Qgis was unable to determine the type and srid of column </source>
        <translation type="unfinished">Qgis nedokázal rozpoznať typ a srid stĺpca </translation>
    </message>
    <message>
        <source>Unable to determine table access privileges for the </source>
        <translation type="unfinished">Nemožno určiť práva prístupu k tabuľke pre </translation>
    </message>
    <message>
        <source>Error while adding features</source>
        <translation>Chyba pri pridávaní objektov</translation>
    </message>
    <message>
        <source>Error while deleting features</source>
        <translation>Chyba pri mazaní objektov</translation>
    </message>
    <message>
        <source>Error while adding attributes</source>
        <translation>Chyba pri pridávaní atribútov</translation>
    </message>
    <message>
        <source>Error while deleting attributes</source>
        <translation>Chyba pri mazaní atribútov</translation>
    </message>
    <message>
        <source>Error while changing attributes</source>
        <translation type="unfinished">Chyba pri úprave atribútov</translation>
    </message>
    <message>
        <source>Error while changing geometry values</source>
        <translation>Chyba pri zmene hodnôt geometrie</translation>
    </message>
    <message>
        <source>unexpected PostgreSQL error</source>
        <translation>neočakávaná chyba PostgreSQL</translation>
    </message>
</context>
<context>
    <name>QgsPostgresProvider::Conn</name>
    <message>
        <source>No GEOS Support!</source>
        <translation>Žiadna podpora GEOSu!</translation>
    </message>
    <message>
        <source>Your PostGIS installation has no GEOS support.
Feature selection and identification will not work properly.
Please install PostGIS with GEOS support (http://geos.refractions.net)</source>
        <translation>Vaša inštalácia PostGIS nobsahuje GEOS.
Výber objektov a identifikácia nebudú pracovať správne.
Prosím nainštalujte PostGIS s podporou GEOSu (http://geos.refractions.net)</translation>
    </message>
</context>
<context>
    <name>QgsProjectPropertiesBase</name>
    <message>
        <source>Project Properties</source>
        <translation>Vlastnosti projektu</translation>
    </message>
    <message>
        <source>Meters</source>
        <translation>Metre</translation>
    </message>
    <message>
        <source>Feet</source>
        <translation>Stopy</translation>
    </message>
    <message>
        <source>Decimal degrees</source>
        <translation>Desiatkové stupne</translation>
    </message>
    <message>
        <source>Default project title</source>
        <translation>Štandardný názov projektu</translation>
    </message>
    <message>
        <source>Automatic</source>
        <translation>Automaticky</translation>
    </message>
    <message>
        <source>Automatically sets the number of decimal places in the mouse position display</source>
        <translation>Automaticky nastaví počet počet desatinných miest v zobrazovaní polohy myši (kurzora)</translation>
    </message>
    <message>
        <source>The number of decimal places that are used when displaying the mouse position is automatically set to be enough so that moving the mouse by one pixel gives a change in the position display</source>
        <translation>Počet desatinných miest, ktoré sú použité pri zobrazovaní polohy myši je automaticky nastavený tak, aby posun myši o jeden pixel spôsobil zmenu v zobrazovanej polohe </translation>
    </message>
    <message>
        <source>Manual</source>
        <translation>Ručne</translation>
    </message>
    <message>
        <source>Sets the number of decimal places to use for the mouse position display</source>
        <translation>Nastaví počet desatinných miest použitých pri zobrazovaní polohy myši</translation>
    </message>
    <message>
        <source>The number of decimal places for the manual option</source>
        <translation>Počet desatinných miest pre ručnú voľbu</translation>
    </message>
    <message>
        <source>decimal places</source>
        <translation>desatinné miesta</translation>
    </message>
    <message>
        <source>Projection</source>
        <translation type="obsolete">Mapové zobrazenie</translation>
    </message>
    <message>
        <source>Enable on the fly projection</source>
        <translation type="obsolete">Povoliť priamy prevod medzi zobrazeniami</translation>
    </message>
    <message>
        <source>General</source>
        <translation>Všeobecné</translation>
    </message>
    <message>
        <source>Digitizing</source>
        <translation>Digitalizácia</translation>
    </message>
    <message>
        <source>Precision</source>
        <translation>Presnosť</translation>
    </message>
    <message>
        <source>Descriptive project name</source>
        <translation>Popisný názov projektu</translation>
    </message>
    <message>
        <source>Enable topological editing</source>
        <translation type="unfinished">Povoliť topologické úpravy</translation>
    </message>
    <message>
        <source>Snapping options...</source>
        <translation>Nastavenie zameriavania...</translation>
    </message>
    <message>
        <source>Avoid intersections of new polygons</source>
        <translation>Vyhnúť sa prieniku nových polygónov</translation>
    </message>
    <message>
        <source>Title and colors</source>
        <translation>Názov a farby</translation>
    </message>
    <message>
        <source>Project title</source>
        <translation>Názov projektu</translation>
    </message>
    <message>
        <source>Selection color</source>
        <translation>Farba výberu</translation>
    </message>
    <message>
        <source>Background color</source>
        <translation>Farba pozadia</translation>
    </message>
    <message>
        <source>Map units</source>
        <translation>Mapové jednotky</translation>
    </message>
    <message>
        <source>Coordinate Reference System (CRS)</source>
        <translation>Referenčný súradnicový systém (CRS)</translation>
    </message>
    <message>
        <source>Enable &apos;on the fly&apos; CRS transformation</source>
        <translation>Povoliť priamy prevod medzi súradnicovými systémami</translation>
    </message>
</context>
<context>
    <name>QgsProjectionSelector</name>
    <message>
        <source>QGIS SRSID: </source>
        <translation type="obsolete">QGIS SRSID:</translation>
    </message>
    <message>
        <source>PostGIS SRID: </source>
        <translation type="obsolete">Postgis SRID:</translation>
    </message>
    <message>
        <source>User Defined Coordinate Systems</source>
        <translation>Súradnicové systémy definované užívateľom</translation>
    </message>
    <message>
        <source>Geographic Coordinate Systems</source>
        <translation>Zemepisné súradnicové systémy</translation>
    </message>
    <message>
        <source>Projected Coordinate Systems</source>
        <translation>Mapové súradnicové systémy</translation>
    </message>
    <message>
        <source>Resource Location Error</source>
        <translation>Chyba pri lokalizácii zdrojov</translation>
    </message>
    <message>
        <source>Error reading database file from: 
 %1
Because of this the projection selector will not work...</source>
        <translation>Chyba pri čítaní databázového súboru: 
 %1
Kvôli tomu nebude fungovať výber mapového zobrazenia...</translation>
    </message>
</context>
<context>
    <name>QgsProjectionSelectorBase</name>
    <message>
        <source>Projection Selector</source>
        <translation type="obsolete">Výber zobrazenia</translation>
    </message>
    <message>
        <source>Projection</source>
        <translation type="obsolete">Zobrazenie</translation>
    </message>
    <message>
        <source>Search</source>
        <translation>Vyhľadávanie</translation>
    </message>
    <message>
        <source>Find</source>
        <translation>Nájsť</translation>
    </message>
    <message>
        <source>Postgis SRID</source>
        <translation type="obsolete">Postgis SRID</translation>
    </message>
    <message>
        <source>EPSG ID</source>
        <translation>EPSG ID</translation>
    </message>
    <message>
        <source>QGIS SRSID</source>
        <translation type="obsolete">QGIS SRSID</translation>
    </message>
    <message>
        <source>Name</source>
        <translation>Názov zobrazenia</translation>
    </message>
    <message>
        <source>Spatial Reference System</source>
        <translation type="obsolete">Priestorový referenčný systém</translation>
    </message>
    <message>
        <source>Id</source>
        <translation type="obsolete">Id</translation>
    </message>
    <message>
        <source>Coordinate Reference System Selector</source>
        <translation>Výber referenčného súradnicového systému</translation>
    </message>
    <message>
        <source>Coordinate Reference System</source>
        <translation>Referenčný súradnicový systém</translation>
    </message>
    <message>
        <source>EPSG</source>
        <translation>EPSG</translation>
    </message>
    <message>
        <source>ID</source>
        <translation>ID</translation>
    </message>
</context>
<context>
    <name>QgsPythonDialog</name>
    <message>
        <source>Python console</source>
        <translation>Konzola Pythonu</translation>
    </message>
    <message>
        <source>&gt;&gt;&gt;</source>
        <translation>&gt;&gt;&gt;</translation>
    </message>
    <message>
        <source>To access Quantum GIS environment from this python console use object from global scope which is an instance of QgisInterface class.&lt;br&gt;Usage e.g.: iface.zoomFull()</source>
        <translation type="unfinished">Na prístup k prostrediu Quantum GIS z tejto konzoly Pythonu použite objekt z global scope, ktorý je inštanciou triedy QgisInterface. &lt;br&gt;Použitie: napr. iface.zoomFull()</translation>
    </message>
</context>
<context>
    <name>QgsQuickPrint</name>
    <message>
        <source> km</source>
        <translation> km</translation>
    </message>
    <message>
        <source> mm</source>
        <translation> mm</translation>
    </message>
    <message>
        <source> cm</source>
        <translation> cm</translation>
    </message>
    <message>
        <source> m</source>
        <translation> m</translation>
    </message>
    <message>
        <source> miles</source>
        <translation> míľ</translation>
    </message>
    <message>
        <source> mile</source>
        <translation> míľa</translation>
    </message>
    <message>
        <source> inches</source>
        <translation> palcov</translation>
    </message>
    <message>
        <source> foot</source>
        <translation> stopa</translation>
    </message>
    <message>
        <source> feet</source>
        <translation> stôp</translation>
    </message>
    <message>
        <source> degree</source>
        <translation> stupeň</translation>
    </message>
    <message>
        <source> degrees</source>
        <translation>stupňov</translation>
    </message>
    <message>
        <source> unknown</source>
        <translation> neznáme jednotky</translation>
    </message>
</context>
<context>
    <name>QgsRasterLayer</name>
    <message>
        <source>Not Set</source>
        <translation>Nenastavené</translation>
    </message>
    <message>
        <source>Driver:</source>
        <translation>Ovládač:</translation>
    </message>
    <message>
        <source>Dimensions:</source>
        <translation>Rozmery:</translation>
    </message>
    <message>
        <source>X: </source>
        <translation type="obsolete">X:</translation>
    </message>
    <message>
        <source> Y: </source>
        <translation type="obsolete">Y:</translation>
    </message>
    <message>
        <source> Bands: </source>
        <translation type="obsolete">Kanály: </translation>
    </message>
    <message>
        <source>Origin:</source>
        <translation>Počiatok:</translation>
    </message>
    <message>
        <source>Pixel Size:</source>
        <translation>Veľkosť pixela:</translation>
    </message>
    <message>
        <source>Raster Extent: </source>
        <translation type="obsolete">Rozsah rastra:</translation>
    </message>
    <message>
        <source>Clipped area: </source>
        <translation type="obsolete"> Vystrihnutá oblasť:</translation>
    </message>
    <message>
        <source>Pyramid overviews:</source>
        <translation>Prehľad pyramíd:</translation>
    </message>
    <message>
        <source>Band</source>
        <translation>Kanál</translation>
    </message>
    <message>
        <source>Band No</source>
        <translation type="unfinished">Kanál č. </translation>
    </message>
    <message>
        <source>No Stats</source>
        <translation>Žiadna štatistika</translation>
    </message>
    <message>
        <source>No stats collected yet</source>
        <translation>Zatiaľ neboli zozbierané žiadne štatistické údaje</translation>
    </message>
    <message>
        <source>Min Val</source>
        <translation>Minimálna hodnota</translation>
    </message>
    <message>
        <source>Max Val</source>
        <translation>Maximálna hodnota</translation>
    </message>
    <message>
        <source>Range</source>
        <translation>Rozpätie (valencia)</translation>
    </message>
    <message>
        <source>Mean</source>
        <translation>Priemer</translation>
    </message>
    <message>
        <source>Sum of squares</source>
        <translation>Rozptyl</translation>
    </message>
    <message>
        <source>Standard Deviation</source>
        <translation>Štandardná odchýlka</translation>
    </message>
    <message>
        <source>Sum of all cells</source>
        <translation>Suma všetkých buniek</translation>
    </message>
    <message>
        <source>Cell Count</source>
        <translation>Počet buniek</translation>
    </message>
    <message>
        <source>Data Type:</source>
        <translation>Typ údajov:</translation>
    </message>
    <message>
        <source>GDT_Byte - Eight bit unsigned integer</source>
        <translation>GDT_Byte - 8 bitový unsigned integer</translation>
    </message>
    <message>
        <source>GDT_UInt16 - Sixteen bit unsigned integer </source>
        <translation>GDT_UInt16 - 16 bitový unsigned integer </translation>
    </message>
    <message>
        <source>GDT_Int16 - Sixteen bit signed integer </source>
        <translation>GDT_Int16 - 16 bitový signed integer</translation>
    </message>
    <message>
        <source>GDT_UInt32 - Thirty two bit unsigned integer </source>
        <translation>GDT_UInt32 - 32 bitový unsigned integer</translation>
    </message>
    <message>
        <source>GDT_Int32 - Thirty two bit signed integer </source>
        <translation>GDT_Int32 - 32 bitový signed integer</translation>
    </message>
    <message>
        <source>GDT_Float32 - Thirty two bit floating point </source>
        <translation>GDT_Float32 - 32 bitový floating point</translation>
    </message>
    <message>
        <source>GDT_Float64 - Sixty four bit floating point </source>
        <translation>GDT_Float64 - 64 bitový floating point</translation>
    </message>
    <message>
        <source>GDT_CInt16 - Complex Int16 </source>
        <translation>GDT_CInt16 - Complex Int16</translation>
    </message>
    <message>
        <source>GDT_CInt32 - Complex Int32 </source>
        <translation>GDT_CInt32 - Complex Int32</translation>
    </message>
    <message>
        <source>GDT_CFloat32 - Complex Float32 </source>
        <translation>GDT_CFloat32 - Complex Float32</translation>
    </message>
    <message>
        <source>GDT_CFloat64 - Complex Float64 </source>
        <translation>GDT_CFloat64 - Complex Float64</translation>
    </message>
    <message>
        <source>Could not determine raster data type.</source>
        <translation>Nemožno rozpoznať typ rastrových údajov.</translation>
    </message>
    <message>
        <source>Average Magphase</source>
        <translation type="unfinished">Priemerná Magphase</translation>
    </message>
    <message>
        <source>Average</source>
        <translation>Priemer</translation>
    </message>
    <message>
        <source>Layer Spatial Reference System: </source>
        <translation>Referenčný priestorový systém vrstvy: </translation>
    </message>
    <message>
        <source>out of extent</source>
        <translation>mimo rozsahu</translation>
    </message>
    <message>
        <source>null (no data)</source>
        <translation>null (žiadne údaje)</translation>
    </message>
    <message>
        <source>Dataset Description</source>
        <translation>Popis súpravy údajov</translation>
    </message>
    <message>
        <source>No Data Value</source>
        <translation>Údaj bez hodnoty</translation>
    </message>
    <message>
        <source>and all other files</source>
        <translation type="obsolete">a všetky ostatné súbory</translation>
    </message>
    <message>
        <source>NoDataValue not set</source>
        <translation>Hodnota pre prázdne bunky nie je nastavená</translation>
    </message>
    <message>
        <source>Band %1</source>
        <translation>Kanál %1</translation>
    </message>
    <message>
        <source>%1 and all other files (*)</source>
        <translation type="unfinished">%1 a všetky ostatné súbory (*)</translation>
    </message>
    <message>
        <source>Band%1</source>
        <translation type="unfinished">Kanál%1</translation>
    </message>
    <message>
        <source>X: %1 Y: %2 Bands: %3</source>
        <translation type="unfinished">X:%1 Y:%2 Kanály: %3</translation>
    </message>
    <message>
        <source>Project Spatial Reference System: </source>
        <translation type="unfinished">Priestorový referenčný systém projektu:</translation>
    </message>
</context>
<context>
    <name>QgsRasterLayerProperties</name>
    <message>
        <source>Grayscale</source>
        <translation>Odtiene šedej</translation>
    </message>
    <message>
        <source>Pseudocolor</source>
        <translation>Nepravé farby</translation>
    </message>
    <message>
        <source>Freak Out</source>
        <translation type="unfinished">Freak Out</translation>
    </message>
    <message>
        <source>Palette</source>
        <translation type="obsolete">Paleta</translation>
    </message>
    <message>
        <source>Not Set</source>
        <translation>Nenastavené</translation>
    </message>
    <message>
        <source>Columns: </source>
        <translation>Stĺpce: </translation>
    </message>
    <message>
        <source>Rows: </source>
        <translation>Riadky: </translation>
    </message>
    <message>
        <source>No-Data Value: </source>
        <translation>Prázdna hodnota: </translation>
    </message>
    <message>
        <source>n/a</source>
        <translation>--</translation>
    </message>
    <message>
        <source>Write access denied</source>
        <translation>Prístup k zápisu odmietnutý</translation>
    </message>
    <message>
        <source>Write access denied. Adjust the file permissions and try again.

</source>
        <translation type="unfinished">Prístup k zápisu bol odmietnutý. Rozšírte práva súboru a skúste znova.

</translation>
    </message>
    <message>
        <source>Building pyramids failed.</source>
        <translation>Tvorba pyramíd zlyhala.</translation>
    </message>
    <message>
        <source>The file was not writeable. Some formats can not be written to, only read. You can also try to check the permissions and then try again.</source>
        <translation type="obsolete">Súbor nie je zapisovateľný. Niektoré formáty nemožno zapisovať, iba čítať. Skontrolujte práva a skúste znova.</translation>
    </message>
    <message>
        <source>Building pyramid overviews is not supported on this type of raster.</source>
        <translation>Pri tomto type rastra nie je možné zostaviť pyramídy.</translation>
    </message>
    <message>
        <source>Custom Colormap</source>
        <translation type="obsolete">Vlastná farebná škála</translation>
    </message>
    <message>
        <source>No Stretch</source>
        <translation type="unfinished">Bez úprav</translation>
    </message>
    <message>
        <source>Stretch To MinMax</source>
        <translation type="unfinished">Natiahnuť na interval MinMax</translation>
    </message>
    <message>
        <source>Stretch And Clip To MinMax</source>
        <translation type="unfinished">Zúžiť rozsah a natiahnuť na MinMax</translation>
    </message>
    <message>
        <source>Clip To MinMax</source>
        <translation type="unfinished">Zúžiť rozsah na MinMax</translation>
    </message>
    <message>
        <source>Discrete</source>
        <translation type="unfinished">Diskrétna</translation>
    </message>
    <message>
        <source>Linearly</source>
        <translation type="obsolete">Lineárna</translation>
    </message>
    <message>
        <source>Equal interval</source>
        <translation type="unfinished">Rovnaký interval</translation>
    </message>
    <message>
        <source>Quantiles</source>
        <translation>Kvantily</translation>
    </message>
    <message>
        <source>Description</source>
        <translation>Popis</translation>
    </message>
    <message>
        <source>Large resolution raster layers can slow navigation in QGIS.</source>
        <translation type="unfinished">Rastrové vrstvy vo vyššom rozlíšení môžu spomaliť navigáciu v QGISe.</translation>
    </message>
    <message>
        <source>By creating lower resolution copies of the data (pyramids) performance can be considerably improved as QGIS selects the most suitable resolution to use depending on the level of zoom.</source>
        <translation type="unfinished">Vytvorením kópií údajov v nížšom rozlíšení (pyramíd) možno badateľne zvýšiť výkon tak, že QGIS pri vykresľovaní rastra vyberie raster s najvhodnejším rozlíšením v závislosti od úrovne pohľadu.</translation>
    </message>
    <message>
        <source>You must have write access in the directory where the original data is stored to build pyramids.</source>
        <translation type="unfinished">Pri tvorbe pyramíd je potrebné mať právo na zápis do adresára, v ktorom sa nachádza originálny súbor.</translation>
    </message>
    <message>
        <source>Please note that building pyramids may alter the original data file and once created they cannot be removed!</source>
        <translation type="obsolete">Pri tvorbe pyramíd treba mať na pamäti, že pyramídy môžu nahradiť pôvodný súbor údajov a pokiaľ sú už raz vytvorené, nemožno ich odobrať!</translation>
    </message>
    <message>
        <source>Please note that building pyramids could corrupt your image - always make a backup of your data first!</source>
        <translation type="obsolete">Pri tvorbe pyramíd treba brať do úvahy, že môžu poškodiť pôvodný súbor, preto si vždy pred touto operáciou  vytvorte záložnú kópiu vašich údajov!</translation>
    </message>
    <message>
        <source>Red</source>
        <translation>Červená</translation>
    </message>
    <message>
        <source>Green</source>
        <translation>Zelená</translation>
    </message>
    <message>
        <source>Blue</source>
        <translation>Modrá</translation>
    </message>
    <message>
        <source>Percent Transparent</source>
        <translation>Priehľadnosť v percentách</translation>
    </message>
    <message>
        <source>Gray</source>
        <translation type="unfinished">Šedá</translation>
    </message>
    <message>
        <source>Indexed Value</source>
        <translation type="unfinished">Indexovaná hodnota</translation>
    </message>
    <message>
        <source>User Defined</source>
        <translation>Užívateľsky definované</translation>
    </message>
    <message>
        <source>No Scaling</source>
        <translation type="obsolete">Bez</translation>
    </message>
    <message>
        <source>No-Data Value: Not Set</source>
        <translation>Prázdna hodnota: Nenastavená</translation>
    </message>
    <message>
        <source>Save file</source>
        <translation>Uložiť súbor</translation>
    </message>
    <message>
        <source>Textfile (*.txt)</source>
        <translation>Textový súbor (*.txt)</translation>
    </message>
    <message>
        <source>QGIS Generated Transparent Pixel Value Export File</source>
        <translation type="unfinished">Súbor s priehľadnými hodnotami pixelov vygenerovaný QGISom</translation>
    </message>
    <message>
        <source>Open file</source>
        <translation>Otvoriť súbor</translation>
    </message>
    <message>
        <source>Import Error</source>
        <translation>Chyba pri importe</translation>
    </message>
    <message>
        <source>The following lines contained errors

</source>
        <translation>Nasledujúce riadky obsahujú chyby

</translation>
    </message>
    <message>
        <source>Read access denied</source>
        <translation type="unfinished">Prístup k čítaniu súboru zamietnutý</translation>
    </message>
    <message>
        <source>Read access denied. Adjust the file permissions and try again.

</source>
        <translation type="unfinished">Prírup čítania zamietnutý. Upravte práva súboru a skúste znova.

</translation>
    </message>
    <message>
        <source>Color Ramp</source>
        <translation type="unfinished">Farebná škála</translation>
    </message>
    <message>
        <source>Please note that building internal pyramids may alter the original data file and once created they cannot be removed!</source>
        <translation>Pri tvorbe pyramíd treba mať na pamäti, že pyramídy môžu nahradiť pôvodný súbor údajov a pokiaľ sú už raz vytvorené, nemožno ich odobrať!</translation>
    </message>
    <message>
        <source>Please note that building internal pyramids could corrupt your image - always make a backup of your data first!</source>
        <translation>Pri tvorbe pyramíd treba brať do úvahy, že môžu poškodiť pôvodný súbor, preto si vždy pred touto operáciou  vytvorte záložnú kópiu vašich údajov!</translation>
    </message>
    <message>
        <source>Default</source>
        <translation type="unfinished">Predvolené</translation>
    </message>
    <message>
        <source>The file was not writeable. Some formats do not support pyramid overviews. Consult the GDAL documentation if in doubt.</source>
        <translation>Do tohoto súboru nemožno zapisovať. Niektoré formáty nepodporujú pyramídové náhľady. V prípade ďalších otázok si prezrite dokumentáciu ku knižnici GDAL.</translation>
    </message>
    <message>
        <source>Default Style</source>
        <translation>Predvolený štýl</translation>
    </message>
    <message>
        <source>QGIS Layer Style File (*.qml)</source>
        <translation>QGIS Súbor štýlu vrstvy (*.qml)</translation>
    </message>
    <message>
        <source>Saved Style</source>
        <translation type="unfinished">Uložený štýl</translation>
    </message>
    <message>
        <source>QGIS</source>
        <translation>QGIS</translation>
    </message>
    <message>
        <source>Unknown style format: </source>
        <translation>Neznámy formát štýlu: </translation>
    </message>
    <message>
        <source>Colormap</source>
        <translation type="unfinished">Mapa farieb</translation>
    </message>
    <message>
        <source>Linear</source>
        <translation type="unfinished">Lineárna</translation>
    </message>
    <message>
        <source>Exact</source>
        <translation type="unfinished">Presne</translation>
    </message>
    <message>
        <source>Custom color map entry</source>
        <translation type="unfinished">Vlastná položka mapy farieb</translation>
    </message>
    <message>
        <source>QGIS Generated Color Map Export File</source>
        <translation type="unfinished">QGIS Export farebnej mapy do súboru</translation>
    </message>
    <message>
        <source>Load Color Map</source>
        <translation type="unfinished">Nahrať mapu farieb</translation>
    </message>
    <message>
        <source>The color map for Band %n failed to load</source>
        <translation type="obsolete">
        </translation>
    </message>
    <message>
        <source>Building internal pyramid overviews is not supported on raster layers with JPEG compression.</source>
        <translation type="unfinished">Tvorba pyramídových náhľadov nie je podporovaná pre rastrové vrstvy s kompresiou JPEG.</translation>
    </message>
    <message>
        <source>Note: Minimum Maximum values are estimates or user defined</source>
        <translation type="unfinished">Poznámka: Hodnoty maxima a minima sú odhadnuté, alebo užívateľom definované</translation>
    </message>
    <message>
        <source>Note: Minimum Maximum values are actual values computed from the band(s)</source>
        <translation type="unfinished">Poznámka: Hodnoty maxima a minima sú skutočné hondoty vypočítané z kanálu (kanálov)</translation>
    </message>
</context>
<context>
    <name>QgsRasterLayerPropertiesBase</name>
    <message>
        <source>Raster Layer Properties</source>
        <translation>Vlastnosti rastrovej vrstvy</translation>
    </message>
    <message>
        <source>General</source>
        <translation>Všeobecné</translation>
    </message>
    <message>
        <source>Legend:</source>
        <translation type="obsolete"> Legenda:</translation>
    </message>
    <message>
        <source>No Data:</source>
        <translation>Žiadne údaje:</translation>
    </message>
    <message>
        <source>Symbology</source>
        <translation>Symbolika</translation>
    </message>
    <message>
        <source>&lt;p align=&quot;right&quot;&gt;Full&lt;/p&gt;</source>
        <translation>&lt;p align=&quot;right&quot;&gt;Úplná&lt;/p&gt;</translation>
    </message>
    <message>
        <source>None</source>
        <translation>Žiadna</translation>
    </message>
    <message>
        <source>Gray</source>
        <translation type="obsolete">Šedá</translation>
    </message>
    <message>
        <source>Metadata</source>
        <translation>Metaúdaje</translation>
    </message>
    <message>
        <source>Pyramids</source>
        <translation>Pyramídy</translation>
    </message>
    <message>
        <source>Average</source>
        <translation>Priemer</translation>
    </message>
    <message>
        <source>Nearest Neighbour</source>
        <translation>Metóda najbližšieho suseda</translation>
    </message>
    <message>
        <source>Thumbnail</source>
        <translation>Miniatúry</translation>
    </message>
    <message>
        <source>Columns:</source>
        <translation>Stĺpce:</translation>
    </message>
    <message>
        <source>Rows:</source>
        <translation>Riadky:</translation>
    </message>
    <message>
        <source>Palette:</source>
        <translation type="obsolete">Paleta:</translation>
    </message>
    <message>
        <source>Maximum scale at which this layer will be displayed. </source>
        <translation>Maximálna mierka pri ktorej bude táto vrstva zobrazená.</translation>
    </message>
    <message>
        <source>Minimum scale at which this layer will be displayed. </source>
        <translation>Minimálna mierka pri ktorej bude táto vrstva zobrazená.</translation>
    </message>
    <message>
        <source>Histogram</source>
        <translation>Stĺpcový graf (histogram)</translation>
    </message>
    <message>
        <source>Options</source>
        <translation>Voľby</translation>
    </message>
    <message>
        <source>Chart Type</source>
        <translation>Typ grafu</translation>
    </message>
    <message>
        <source>Refresh</source>
        <translation>Obnoviť</translation>
    </message>
    <message>
        <source>Change</source>
        <translation type="obsolete">Zmeniť</translation>
    </message>
    <message>
        <source>Max</source>
        <translation type="unfinished">Max</translation>
    </message>
    <message>
        <source>Min</source>
        <translation type="unfinished">Min</translation>
    </message>
    <message>
        <source>&lt;b&gt;&lt;font color=&apos;blue&apos;&gt;Max&lt;/font&gt;&lt;/b&gt;</source>
        <translation type="obsolete">&lt;b&gt;&lt;font color=&apos;blue&apos;&gt;Max&lt;/font&gt;&lt;/b&gt;</translation>
    </message>
    <message>
        <source>&lt;b&gt;&lt;font color=&apos;green&apos;&gt;Min&lt;/font&gt;&lt;/b&gt;</source>
        <translation type="obsolete">&lt;b&gt;&lt;font color=&apos;green&apos;&gt;Min&lt;/font&gt;&lt;/b&gt;</translation>
    </message>
    <message>
        <source>&lt;b&gt;&lt;font color=&apos;green&apos;&gt;Max&lt;/font&gt;&lt;/b&gt;</source>
        <translation type="obsolete">&lt;b&gt;&lt;font color=&apos;green&apos;&gt;Max&lt;/font&gt;&lt;/b&gt;</translation>
    </message>
    <message>
        <source>&lt;b&gt;&lt;font color=&apos;red&apos;&gt;Min&lt;/font&gt;&lt;/b&gt;</source>
        <translation type="obsolete">&lt;b&gt;&lt;font color=&apos;red&apos;&gt;Min&lt;/font&gt;&lt;/b&gt;</translation>
    </message>
    <message>
        <source>&lt;b&gt;&lt;font color=&apos;red&apos;&gt;Max&lt;/font&gt;&lt;/b&gt;</source>
        <translation type="obsolete">&lt;b&gt;&lt;font color=&apos;red&apos;&gt;Max&lt;/font&gt;&lt;/b&gt;</translation>
    </message>
    <message>
        <source>&lt;b&gt;&lt;font color=&apos;blue&apos;&gt;Min&lt;/font&gt;&lt;/b&gt;</source>
        <translation type="obsolete">&lt;b&gt;&lt;font color=&apos;blue&apos;&gt;Min&lt;/font&gt;&lt;/b&gt;</translation>
    </message>
    <message>
        <source>&lt;b&gt;&lt;font color=&apos;blue&apos;&gt;Blue&lt;/font&gt;&lt;/b&gt;</source>
        <translation type="obsolete">&lt;b&gt;&lt;font color=&apos;blue&apos;&gt;Modrá&lt;/font&gt;&lt;/b&gt;</translation>
    </message>
    <message>
        <source>&lt;b&gt;&lt;font color=&apos;green&apos;&gt;Green&lt;/font&gt;&lt;/b&gt;</source>
        <translation type="obsolete">&lt;b&gt;&lt;font color=&apos;green&apos;&gt;Zelená&lt;/font&gt;&lt;/b&gt;</translation>
    </message>
    <message>
        <source>&lt;b&gt;&lt;font color=&apos;red&apos;&gt;Red&lt;/font&gt;&lt;/b&gt;</source>
        <translation type="obsolete">&lt;b&gt;&lt;font color=&apos;red&apos;&gt;Červená&lt;/font&gt;&lt;/b&gt;</translation>
    </message>
    <message>
        <source> 00%</source>
        <translation type="unfinished"> 00%</translation>
    </message>
    <message>
        <source>Render as</source>
        <translation type="unfinished">Vykreslovať ako</translation>
    </message>
    <message>
        <source>...</source>
        <translation>...</translation>
    </message>
    <message>
        <source>Colormap</source>
        <translation type="unfinished">Farebná škála</translation>
    </message>
    <message>
        <source>Delete entry</source>
        <translation>Zmazať položku</translation>
    </message>
    <message>
        <source>Classify</source>
        <translation type="unfinished">Klasifikovať</translation>
    </message>
    <message>
        <source>1</source>
        <translation type="unfinished">1</translation>
    </message>
    <message>
        <source>2</source>
        <translation type="unfinished">2</translation>
    </message>
    <message>
        <source>Three band color</source>
        <translation type="unfinished">Trojknálová farebná</translation>
    </message>
    <message>
        <source>Single band gray</source>
        <translation type="unfinished">Jednokanálová v odtieňoch šedi</translation>
    </message>
    <message>
        <source>Color map</source>
        <translation type="unfinished">Mapa farieb</translation>
    </message>
    <message>
        <source>Invert color map</source>
        <translation type="unfinished">Invertovať mapu farieb</translation>
    </message>
    <message>
        <source>RGB mode band selection</source>
        <translation type="obsolete">RGB režim výberu kanála</translation>
    </message>
    <message>
        <source>Grayscale band selection</source>
        <translation type="obsolete">Výber odtieňov šedi</translation>
    </message>
    <message>
        <source>RGB scaling</source>
        <translation type="obsolete">Škálovanie RGB</translation>
    </message>
    <message>
        <source>Std. deviation</source>
        <translation type="obsolete">Štd. odchýlka</translation>
    </message>
    <message>
        <source>Custom min / max values</source>
        <translation type="unfinished">Vlastné min / max hodnoty</translation>
    </message>
    <message>
        <source>Grayscale band scaling</source>
        <translation type="obsolete">Škálovanie odtieňov šedi</translation>
    </message>
    <message>
        <source>Load min / max values from band</source>
        <translation type="unfinished">Nahrať z kanála min / max hodnotu</translation>
    </message>
    <message>
        <source>Estimate (faster)</source>
        <translation type="unfinished">Odhad (rýchlejšie)</translation>
    </message>
    <message>
        <source>Load</source>
        <translation type="unfinished">Nahrať</translation>
    </message>
    <message>
        <source>Actual (slower)</source>
        <translation type="unfinished">Presné (pomalšie)</translation>
    </message>
    <message>
        <source>Contrast enhancement</source>
        <translation type="unfinished">Pridanie kontrastu</translation>
    </message>
    <message>
        <source>Current</source>
        <translation type="unfinished">Aktuálne</translation>
    </message>
    <message>
        <source>Save current contrast enhancement algorithm as default. This setting will be persistent between QGIS sessions.</source>
        <translation type="unfinished">Uloží práve nastavený algoritmus zväčšenia kontrastu ako predvolený. Toto nastavenie sa zachová medzi sedeniami QGISu.</translation>
    </message>
    <message>
        <source>Saves current contrast enhancement algorithm as a default. This setting will be persistent between QGIS sessions.</source>
        <translation type="unfinished">Uloží práve nastavený algoritmus zväčšenia kontrastu ako predvolený. Toto nastavenie sa zachová medzi sedeniami QGISu.</translation>
    </message>
    <message>
        <source>Default</source>
        <translation type="unfinished">Predvolené</translation>
    </message>
    <message>
        <source>TextLabel</source>
        <translation type="unfinished">TextLabel</translation>
    </message>
    <message>
        <source>Transparency</source>
        <translation type="unfinished">Priehľadnosť</translation>
    </message>
    <message>
        <source>Global transparency</source>
        <translation type="unfinished">Globálna priehľadnosť</translation>
    </message>
    <message>
        <source>No data value</source>
        <translation>Hodnota prázdnej bunky</translation>
    </message>
    <message>
        <source>Reset no data value</source>
        <translation type="unfinished">Znovunastaviť hodnotu prázdnej bunky</translation>
    </message>
    <message>
        <source>Custom transparency options</source>
        <translation type="unfinished">Vlastné nastavenia priehľadnosti</translation>
    </message>
    <message>
        <source>Transparency band</source>
        <translation type="unfinished">Kanál priehľadnosti</translation>
    </message>
    <message>
        <source>Transparency layer;</source>
        <translation type="obsolete">Vrstva priehľadnosti;</translation>
    </message>
    <message>
        <source>Transparent pixel list</source>
        <translation>Zoznam priehľadných pixelov</translation>
    </message>
    <message>
        <source>Add values manually</source>
        <translation type="unfinished">Pridať hodnoty ručne</translation>
    </message>
    <message>
        <source>Add Values from display</source>
        <translation type="unfinished">Pridať hodnoty z obrazovky</translation>
    </message>
    <message>
        <source>Remove selected row</source>
        <translation type="unfinished">Odobrať vybraný riadok</translation>
    </message>
    <message>
        <source>Default values</source>
        <translation type="unfinished">Predvolené hodnoty</translation>
    </message>
    <message>
        <source>Import from file</source>
        <translation type="unfinished">Import zo súboru</translation>
    </message>
    <message>
        <source>Export to file</source>
        <translation type="unfinished">Export do súboru</translation>
    </message>
    <message>
        <source>Number of entries</source>
        <translation type="unfinished">Počet položiek</translation>
    </message>
    <message>
        <source>Color interpolation</source>
        <translation type="unfinished">Interpolácia farieb</translation>
    </message>
    <message>
        <source>Classification mode</source>
        <translation type="unfinished">Režim klasifikácie</translation>
    </message>
    <message>
        <source>Spatial reference system</source>
        <translation type="obsolete">Priestorový referenčný systém</translation>
    </message>
    <message>
        <source>Scale dependent visibility</source>
        <translation>Viditeľnosť závislá od mierky</translation>
    </message>
    <message>
        <source>Maximum</source>
        <translation type="unfinished">Maximum</translation>
    </message>
    <message>
        <source>Minimum</source>
        <translation type="unfinished">Minimum</translation>
    </message>
    <message>
        <source>Show debug info</source>
        <translation type="obsolete">Zobraziť ladiace informácie</translation>
    </message>
    <message>
        <source>Layer source</source>
        <translation>Zdroj vrstvy</translation>
    </message>
    <message>
        <source>Display name</source>
        <translation type="unfinished">Zobraziť meno</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Lucida Grande&apos;; font-size:13pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot;-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Lucida Grande&apos;; font-size:13pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot;-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>Pyramid resolutions</source>
        <translation>Rozlíšenie pyramíd</translation>
    </message>
    <message>
        <source>Resampling method</source>
        <translation>Metóda prevzorkovania</translation>
    </message>
    <message>
        <source>Build pyramids</source>
        <translation>Vytvoriť pyramídy</translation>
    </message>
    <message>
        <source>Line graph</source>
        <translation>Čiarový graf</translation>
    </message>
    <message>
        <source>Bar chart</source>
        <translation>Stĺpcový graf</translation>
    </message>
    <message>
        <source>Column count</source>
        <translation>Počet stĺpcov</translation>
    </message>
    <message>
        <source>Out of range OK?</source>
        <translation>Môže byť mimo rozsahu?</translation>
    </message>
    <message>
        <source>Allow approximation</source>
        <translation>Povoliť odhad</translation>
    </message>
    <message>
        <source>RGB mode band selection and scaling</source>
        <translation type="unfinished">RGB režim výberu a škálovania</translation>
    </message>
    <message>
        <source>Red band</source>
        <translation>Červený kanál</translation>
    </message>
    <message>
        <source>Green band</source>
        <translation>Zelený kanál</translation>
    </message>
    <message>
        <source>Blue band</source>
        <translation>Modrý kanál</translation>
    </message>
    <message>
        <source>Default R:1 G:2 B:3</source>
        <translation type="unfinished">Predvolene: R:1 G:2 B:3</translation>
    </message>
    <message>
        <source>Save current band combination as default. This setting will be persistent between QGIS sessions.</source>
        <translation type="obsolete">Uložiť aktuálnu kombináciu kanálov ako predvolenú. Toto nastavenia ostane zaochované pre ďalšie sedenia QGISu.</translation>
    </message>
    <message>
        <source>Save current band combination as a default. This setting will be persistent between QGIS sessions.</source>
        <translation type="obsolete">Uloží aktuálnu kombináciu kanálov ako predvolenú. Toto nastavenie ostane zachované medzi sedeniami QGISu.</translation>
    </message>
    <message>
        <source>Red min</source>
        <translation type="unfinished">Červená min</translation>
    </message>
    <message>
        <source>Red max</source>
        <translation type="unfinished">Červená max</translation>
    </message>
    <message>
        <source>Green min</source>
        <translation type="unfinished">Zelená min</translation>
    </message>
    <message>
        <source>Green max</source>
        <translation type="unfinished">Zelená max</translation>
    </message>
    <message>
        <source>Blue min</source>
        <translation type="unfinished">Modrá min</translation>
    </message>
    <message>
        <source>Blue max</source>
        <translation type="unfinished">Modrá max</translation>
    </message>
    <message>
        <source>Use standard deviation</source>
        <translation type="unfinished">Použiť štandardnú odchýlku</translation>
    </message>
    <message>
        <source>Single band properties</source>
        <translation type="unfinished">Vlastnosti jednokanála</translation>
    </message>
    <message>
        <source>Gray band</source>
        <translation type="unfinished">Kanál šedej</translation>
    </message>
    <message>
        <source>Note:</source>
        <translation type="unfinished">Poznámka:</translation>
    </message>
    <message>
        <source>Notes</source>
        <translation type="unfinished">Poznámky</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot;-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="obsolete">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Verdana&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot;-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;; font-size:9pt;&quot;&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>Build pyramids internally if possible</source>
        <translation type="unfinished">Pokiaľ je to možné, vytvoriť pyramídy interne</translation>
    </message>
    <message>
        <source>Restore Default Style</source>
        <translation type="unfinished">Obnoviť predvolený štýl</translation>
    </message>
    <message>
        <source>Save As Default</source>
        <translation type="unfinished">Uložiť ako predvolené</translation>
    </message>
    <message>
        <source>Load Style ...</source>
        <translation type="unfinished">Nahrať štýl ...</translation>
    </message>
    <message>
        <source>Save Style ...</source>
        <translation>Uložiť štýl ...</translation>
    </message>
    <message>
        <source>Add entry</source>
        <translation>Pridať položku</translation>
    </message>
    <message>
        <source>Sort</source>
        <translation>Usporiadať</translation>
    </message>
    <message>
        <source>Load color map from band</source>
        <translation type="unfinished">Nahrať farebnú mapu z kanála</translation>
    </message>
    <message>
        <source>Load color map from file</source>
        <translation type="unfinished">Nahrať farebnú mapu zo súboru</translation>
    </message>
    <message>
        <source>Export color map to file</source>
        <translation type="unfinished">Exportovať farebnú mapu do súboru</translation>
    </message>
    <message>
        <source>Generate new color map</source>
        <translation type="unfinished">Vygenerovať novú mapu farieb</translation>
    </message>
    <message>
        <source>Coordinate reference system</source>
        <translation>Referenčný súradnicový systém</translation>
    </message>
    <message>
        <source>Change ...</source>
        <translation>Zmeniť ...</translation>
    </message>
    <message>
        <source>Legend</source>
        <translation>Legenda</translation>
    </message>
    <message>
        <source>Palette</source>
        <translation>Paleta</translation>
    </message>
    <message>
        <source>&lt;!DOCTYPE HTML PUBLIC &quot;-//W3C//DTD HTML 4.0//EN&quot; &quot;http://www.w3.org/TR/REC-html40/strict.dtd&quot;&gt;
&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;DejaVu Sans&apos;; font-size:9pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot;-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;;&quot;&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="unfinished">&lt;!DOCTYPE HTML PUBLIC &quot;-//W3C//DTD HTML 4.0//EN&quot; &quot;http://www.w3.org/TR/REC-html40/strict.dtd&quot;&gt;&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;p, li { white-space: pre-wrap; }&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;DejaVu Sans&apos;; font-size:9pt; font-weight:400; font-style:normal;&quot;&gt;&lt;p style=&quot;-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;Sans Serif&apos;;&quot;&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
</context>
<context>
    <name>QgsRunProcess</name>
    <message>
        <source>Unable to run command</source>
        <translation>Nie je možné spustiť príkaz</translation>
    </message>
    <message>
        <source>Starting</source>
        <translation>Štartuje</translation>
    </message>
    <message>
        <source>Done</source>
        <translation>Dokončené</translation>
    </message>
    <message>
        <source>Action</source>
        <translation>Akcia</translation>
    </message>
</context>
<context>
    <name>QgsScaleBarPlugin</name>
    <message>
        <source> metres/km</source>
        <translation> metre/kilometer</translation>
    </message>
    <message>
        <source> feet</source>
        <translation>stopy</translation>
    </message>
    <message>
        <source> degrees</source>
        <translation>stupne</translation>
    </message>
    <message>
        <source> km</source>
        <translation>km</translation>
    </message>
    <message>
        <source> mm</source>
        <translation> mm</translation>
    </message>
    <message>
        <source> cm</source>
        <translation> cm</translation>
    </message>
    <message>
        <source> m</source>
        <translation> m</translation>
    </message>
    <message>
        <source> foot</source>
        <translation> stopa</translation>
    </message>
    <message>
        <source> degree</source>
        <translation> stupeň</translation>
    </message>
    <message>
        <source> unknown</source>
        <translation> neznáma</translation>
    </message>
    <message>
        <source>Top Left</source>
        <translation>Vľavo hore</translation>
    </message>
    <message>
        <source>Bottom Left</source>
        <translation>Vľavo dole</translation>
    </message>
    <message>
        <source>Top Right</source>
        <translation>Vpravo hore</translation>
    </message>
    <message>
        <source>Bottom Right</source>
        <translation>Vpravo dole</translation>
    </message>
    <message>
        <source>Tick Down</source>
        <translation>Čiarka dole</translation>
    </message>
    <message>
        <source>Tick Up</source>
        <translation>Čiarka hore</translation>
    </message>
    <message>
        <source>Bar</source>
        <translation>Lišta</translation>
    </message>
    <message>
        <source>Box</source>
        <translation>Rám</translation>
    </message>
    <message>
        <source>&amp;Decorations</source>
        <translation>&amp;Doplnky</translation>
    </message>
    <message>
        <source>Creates a scale bar that is displayed on the map canvas</source>
        <translation>Vytvorí grafickú mierku zobrazovanú na mapovom plátne</translation>
    </message>
    <message>
        <source>&amp;Scale Bar</source>
        <translation>&amp;Grafická mierka</translation>
    </message>
    <message>
        <source> feet/miles</source>
        <translation type="unfinished"> stopy/míle</translation>
    </message>
    <message>
        <source> miles</source>
        <translation type="unfinished"> míle</translation>
    </message>
    <message>
        <source> mile</source>
        <translation type="unfinished"> míľa</translation>
    </message>
    <message>
        <source> inches</source>
        <translation type="unfinished"> palce</translation>
    </message>
</context>
<context>
    <name>QgsScaleBarPluginGuiBase</name>
    <message>
        <source>Scale Bar Plugin</source>
        <translation>Grafická mierka</translation>
    </message>
    <message>
        <source>Top Left</source>
        <translation>Vľavo hore</translation>
    </message>
    <message>
        <source>Top Right</source>
        <translation>Vpravo hore</translation>
    </message>
    <message>
        <source>Bottom Left</source>
        <translation>Vľavo dole</translation>
    </message>
    <message>
        <source>Bottom Right</source>
        <translation>Vpravo dole</translation>
    </message>
    <message>
        <source>Size of bar:</source>
        <translation>Veľkosť grafickej mierky:</translation>
    </message>
    <message>
        <source>Placement:</source>
        <translation>Umiestnenie:</translation>
    </message>
    <message>
        <source>Tick Down</source>
        <translation>Čiarka dole</translation>
    </message>
    <message>
        <source>Tick Up</source>
        <translation>Čiarka hore</translation>
    </message>
    <message>
        <source>Box</source>
        <translation>Rám</translation>
    </message>
    <message>
        <source>Bar</source>
        <translation>Lišta</translation>
    </message>
    <message>
        <source>Select the style of the scale bar</source>
        <translation>Vybrať vzhľad grafickej mierky</translation>
    </message>
    <message>
        <source>Colour of bar:</source>
        <translation>Farba grafickej mierky:</translation>
    </message>
    <message>
        <source>Scale bar style:</source>
        <translation>Vzhľad grafickej mierky:</translation>
    </message>
    <message>
        <source>Enable scale bar</source>
        <translation>Zapnúť grafickú mierku</translation>
    </message>
    <message>
        <source>Automatically snap to round number on resize</source>
        <translation>Pri zmene veľkosti automaticky skočiť na okrúhle číslo</translation>
    </message>
    <message>
        <source>Click to select the colour</source>
        <translation>Kliknutím vyberiete farbu</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;/head&gt;&lt;body style=&quot; white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;This plugin draws a scale bar on the map. Please note the size option below is a &apos;preferred&apos; size and may have to be altered by QGIS depending on the level of zoom.  The size is measured according to the map units specified in the project properties.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation>  &lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;/head&gt;&lt;body style=&quot; white-space: pre-wrap; font-family:Sans Serif; font-size:9pt; font-weight:400; font-style:normal; text-decoration:none;&quot;&gt;&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Tento zásuvný modul vykreslí na mape grafickú mierku. Prosím nezabúdajte, že tu nastavovaná veľkosť je veľkosť uprednostňovaná, no môže byť nahradená inou v závislosti od nastavenia úrovne pohľadu. Veľkosť je posudzovaná vzhľadom na mapové jednotky nastavené vo Vlastnostiach projektu.&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
</context>
<context>
    <name>QgsSearchQueryBuilder</name>
    <message>
        <source>Found %d matching features.</source>
        <translation type="obsolete">
        </translation>
    </message>
    <message>
        <source>No matching features found.</source>
        <translation>Nenašli sa žiadne zodpovedajúce objekty.</translation>
    </message>
    <message>
        <source>Search results</source>
        <translation>Výsledky vyhľadávania</translation>
    </message>
    <message>
        <source>Search string parsing error</source>
        <translation>Chyba pri kontrole reťazca</translation>
    </message>
    <message>
        <source>No Records</source>
        <translation>Žiadne záznamy</translation>
    </message>
    <message>
        <source>The query you specified results in zero records being returned.</source>
        <translation>Dopyt ktorý ste zadali vrátil nula záznamov.</translation>
    </message>
    <message>
        <source>Search query builder</source>
        <translation type="unfinished">Tvorba dopytov pre vyhľadávanie</translation>
    </message>
    <message>
        <source>Found %1 matching features.</source>
        <translation type="obsolete">Našlo sa %1 zodpovedájúcich objektov.
        </translation>
    </message>
</context>
<context>
    <name>QgsServerSourceSelect</name>
    <message>
        <source>Are you sure you want to remove the </source>
        <translation>Ste si istý, že chcete odstrániť </translation>
    </message>
    <message>
        <source> connection and all associated settings?</source>
        <translation> spojenie a všetky pričlenené nastavenia?</translation>
    </message>
    <message>
        <source>Confirm Delete</source>
        <translation>Potvrdenie mazania</translation>
    </message>
    <message>
        <source>Select Layer</source>
        <translation>Vybrať vrstvu</translation>
    </message>
    <message>
        <source>You must select at least one layer first.</source>
        <translation>Je treba najskôr vybrať aspoň jednu vrstvu.</translation>
    </message>
    <message>
        <source>WMS Provider</source>
        <translation>Nástroj na prístup k údajom WMS</translation>
    </message>
    <message>
        <source>Could not open the WMS Provider</source>
        <translation>Nemožno otvoriť Nástroj na prístup k údajom WMS</translation>
    </message>
    <message>
        <source>Coordinate Reference System (%1 available)</source>
        <translation type="obsolete">Referenčný súradnicový systém (%1 dostupných)
        </translation>
    </message>
    <message>
        <source>Could not understand the response.  The</source>
        <translation type="unfinished">Odpoveď nemožno vyhodnotiť.  </translation>
    </message>
    <message>
        <source>provider said</source>
        <translation type="unfinished">Odpoveď providera</translation>
    </message>
    <message>
        <source>WMS proxies</source>
        <translation type="unfinished">WMS proxy servery</translation>
    </message>
    <message>
        <source>&lt;p&gt;Several WMS servers have been added to the server list. Note that the proxy fields have been left blank and if you access the internet via a web proxy, you will need to individually set the proxy fields with appropriate values.&lt;/p&gt;</source>
        <translation type="obsolete">&lt;p&gt;Niekoľko WMS serverov bolo pridaných do zoznamu serverov. Nezabudnite, že polia pre proxy ostali prázdne a pokiaľ sa pripájate do internetu prostredníctvom proxy servera budete musieť polia proxy s príslušnými hodnotami nastaviť zvlášť.&lt;/p&gt;</translation>
    </message>
    <message>
        <source>Coordinate Reference System</source>
        <translation>Referenčný súradnicový systém</translation>
    </message>
    <message>
        <source>There are no available coordinate reference system for the set of layers you&apos;ve selected.</source>
        <translation type="unfinished">Nie sú dostupné žiadne referenčné súradnicové systémy pre túto množinu vrstiev, ktorú ste vybrali.</translation>
    </message>
    <message>
        <source>Several WMS servers have been added to the server list. Note that if you access the internet via a web proxy, you will need to set the proxy settings in the QGIS options dialog.</source>
        <translation type="unfinished">Niekoľko WMS serverov bolo pridaných do zoznamu serverov. Nezabudnite, že polia pre proxy ostali prázdne a pokiaľ sa pripájate do internetu prostredníctvom proxy servera budete musieť polia proxy s príslušnými hodnotami nastaviť zvlášť.</translation>
    </message>
</context>
<context>
    <name>QgsServerSourceSelectBase</name>
    <message>
        <source>Edit</source>
        <translation>Upraviť</translation>
    </message>
    <message>
        <source>Delete</source>
        <translation>Vymazať</translation>
    </message>
    <message>
        <source>Help</source>
        <translation>Pomocník</translation>
    </message>
    <message>
        <source>F1</source>
        <translation>F1</translation>
    </message>
    <message>
        <source>Add Layer(s) from a Server</source>
        <translation>Pridať vrstvu (vrstvy) zo servera</translation>
    </message>
    <message>
        <source>Server Connections</source>
        <translation>Spojenia</translation>
    </message>
    <message>
        <source>&amp;New</source>
        <translation>&amp;Nové</translation>
    </message>
    <message>
        <source>C&amp;onnect</source>
        <translation>&amp;Pripojiť</translation>
    </message>
    <message>
        <source>Image encoding</source>
        <translation>Kódovanie obrázka</translation>
    </message>
    <message>
        <source>Layers</source>
        <translation>Vrstvy</translation>
    </message>
    <message>
        <source>ID</source>
        <translation>ID</translation>
    </message>
    <message>
        <source>Name</source>
        <translation>Meno</translation>
    </message>
    <message>
        <source>Title</source>
        <translation>Titulok</translation>
    </message>
    <message>
        <source>Abstract</source>
        <translation>Abstrakt</translation>
    </message>
    <message>
        <source>&amp;Add</source>
        <translation>Prid&amp;ať</translation>
    </message>
    <message>
        <source>Alt+A</source>
        <translation>Alt+A</translation>
    </message>
    <message>
        <source>C&amp;lose</source>
        <translation>&amp;Zatvoriť</translation>
    </message>
    <message>
        <source>Alt+L</source>
        <translation>Alt+Z</translation>
    </message>
    <message>
        <source>Ready</source>
        <translation>Pripravený</translation>
    </message>
    <message>
        <source>Coordinate Reference System</source>
        <translation>Referenčný súradnicový systém</translation>
    </message>
    <message>
        <source>Change ...</source>
        <translation>Zmeniť ...</translation>
    </message>
    <message>
        <source>Adds a few example WMS servers</source>
        <translation>Pridať na ukážku niekoľko WMS serverov</translation>
    </message>
    <message>
        <source>Add default servers</source>
        <translation>Pridať prednastavený zoznam serverov</translation>
    </message>
</context>
<context>
    <name>QgsShapeFile</name>
    <message>
        <source>The database gave an error while executing this SQL:</source>
        <translation>Databáza vrátila chyby pri vykonávaní SQL dopytu:</translation>
    </message>
    <message>
        <source>The error was:</source>
        <translation>Chybové hlásenie:</translation>
    </message>
    <message>
        <source>... (rest of SQL trimmed)</source>
        <comment>is appended to a truncated SQL statement</comment>
        <translation>... (zvyšok SQL odrezaný)</translation>
    </message>
    <message>
        <source>Scanning </source>
        <translation type="unfinished">Skenuje sa </translation>
    </message>
</context>
<context>
    <name>QgsSingleSymbolDialog</name>
    <message>
        <source>Solid Line</source>
        <translation>Plná čiara</translation>
    </message>
    <message>
        <source>Dash Line</source>
        <translation>Čiarkovaná čiara</translation>
    </message>
    <message>
        <source>Dot Line</source>
        <translation>Bodkovaná čiara</translation>
    </message>
    <message>
        <source>Dash Dot Line</source>
        <translation>Bodkočiarkovaná čiara</translation>
    </message>
    <message>
        <source>Dash Dot Dot Line</source>
        <translation>Čiarkobodkobodkovaná čiara</translation>
    </message>
    <message>
        <source>No Pen</source>
        <translation>Prázdna čiara</translation>
    </message>
    <message>
        <source>No Brush</source>
        <translation>Bez výplne</translation>
    </message>
    <message>
        <source>Solid</source>
        <translation>Plná</translation>
    </message>
    <message>
        <source>Horizontal</source>
        <translation>Vodorovne</translation>
    </message>
    <message>
        <source>Vertical</source>
        <translation>Zvisle</translation>
    </message>
    <message>
        <source>Cross</source>
        <translation>Krížom</translation>
    </message>
    <message>
        <source>BDiagonal</source>
        <translation>Ulohopriečne zdola</translation>
    </message>
    <message>
        <source>FDiagonal</source>
        <translation>Uhlopriečne zhora</translation>
    </message>
    <message>
        <source>Diagonal X</source>
        <translation>Uhlopriečne X</translation>
    </message>
    <message>
        <source>Dense1</source>
        <translation>Hustá1</translation>
    </message>
    <message>
        <source>Dense2</source>
        <translation>Hustá2</translation>
    </message>
    <message>
        <source>Dense3</source>
        <translation>Hustá3</translation>
    </message>
    <message>
        <source>Dense4</source>
        <translation>Hustá4</translation>
    </message>
    <message>
        <source>Dense5</source>
        <translation>Hustá5</translation>
    </message>
    <message>
        <source>Dense6</source>
        <translation>Hustá6</translation>
    </message>
    <message>
        <source>Dense7</source>
        <translation>Hustá7</translation>
    </message>
    <message>
        <source>Texture</source>
        <translation>Textúra</translation>
    </message>
</context>
<context>
    <name>QgsSingleSymbolDialogBase</name>
    <message>
        <source>Single Symbol</source>
        <translation>Jeden symbol</translation>
    </message>
    <message>
        <source>Size</source>
        <translation>Veľkosť</translation>
    </message>
    <message>
        <source>...</source>
        <translation>...</translation>
    </message>
    <message>
        <source>Point Symbol</source>
        <translation>Bodový symbol</translation>
    </message>
    <message>
        <source>Area scale field</source>
        <translation type="unfinished">Pole s mierkou plochy</translation>
    </message>
    <message>
        <source>Rotation field</source>
        <translation>Pole s rotáciou</translation>
    </message>
    <message>
        <source>Style Options</source>
        <translation>Možnosti štýlu</translation>
    </message>
    <message>
        <source>Outline style</source>
        <translation>Štýl obrysu</translation>
    </message>
    <message>
        <source>Outline color</source>
        <translation>Farba obrysu</translation>
    </message>
    <message>
        <source>Outline width</source>
        <translation>Hrúbka obrysu</translation>
    </message>
    <message>
        <source>Fill color</source>
        <translation>Farba výplne</translation>
    </message>
    <message>
        <source>Fill style</source>
        <translation>Štýl výplne</translation>
    </message>
    <message>
        <source>Label</source>
        <translation>Popis</translation>
    </message>
</context>
<context>
    <name>QgsSnappingDialog</name>
    <message>
        <source>to vertex</source>
        <translation>na uzol</translation>
    </message>
    <message>
        <source>to segment</source>
        <translation>na úsek</translation>
    </message>
    <message>
        <source>to vertex and segment</source>
        <translation>na uzol a úsek</translation>
    </message>
</context>
<context>
    <name>QgsSnappingDialogBase</name>
    <message>
        <source>Snapping options</source>
        <translation>Možnosti zameriavania</translation>
    </message>
    <message>
        <source>Layer</source>
        <translation>Vrstva</translation>
    </message>
    <message>
        <source>Mode</source>
        <translation>Mód</translation>
    </message>
    <message>
        <source>Tolerance</source>
        <translation>Tolerancia</translation>
    </message>
</context>
<context>
    <name>QgsSpit</name>
    <message>
        <source>Are you sure you want to remove the [</source>
        <translation>Naozaj chcete odstrániť spojenie [</translation>
    </message>
    <message>
        <source>] connection and all associated settings?</source>
        <translation>] a všetky s ním previazané nastavenia?</translation>
    </message>
    <message>
        <source>Confirm Delete</source>
        <translation>Potvrdenie mazania</translation>
    </message>
    <message>
        <source>The following Shapefile(s) could not be loaded:

</source>
        <translation>Nasledujúci/e súbor(y) shape nie je možné nahrať:

</translation>
    </message>
    <message>
        <source>REASON: File cannot be opened</source>
        <translation>DÔVOD: Súbor sa nedá otvoriť</translation>
    </message>
    <message>
        <source>REASON: One or both of the Shapefile files (*.dbf, *.shx) missing</source>
        <translation>DÔVOD: Jeden alebo oba súbory shape (*.dbf, *.shx) chýbajú</translation>
    </message>
    <message>
        <source>General Interface Help:</source>
        <translation>Pomocník k hlavnému rozhraniu:</translation>
    </message>
    <message>
        <source>PostgreSQL Connections:</source>
        <translation>Spojenia PostgreSQL:</translation>
    </message>
    <message>
        <source>[New ...] - create a new connection</source>
        <translation>[Nové ...] - vytvoriť nové spojenie</translation>
    </message>
    <message>
        <source>[Edit ...] - edit the currently selected connection</source>
        <translation>[Upraviť ...] - upraviť práve vybraté spojenie</translation>
    </message>
    <message>
        <source>[Remove] - remove the currently selected connection</source>
        <translation>[Odstrániť] - odstrániť práve vybraté spojenie</translation>
    </message>
    <message>
        <source>-you need to select a connection that works (connects properly) in order to import files</source>
        <translation>-aby bol umožnený import súborov je potrebné vybrať fungujúce spojenie (pripojiť sa správne)</translation>
    </message>
    <message>
        <source>-when changing connections Global Schema also changes accordingly</source>
        <translation>-keď sa menia spojenia Global Schema sa mení vzhľadom na ne</translation>
    </message>
    <message>
        <source>Shapefile List:</source>
        <translation>Zoznam súborov shape:</translation>
    </message>
    <message>
        <source>[Add ...] - open a File dialog and browse to the desired file(s) to import</source>
        <translation>[Pridať ...] - otvorí Dialóg so súbormi, kde je možné nalistovať a vybrať požadovaný/é súbor(y) na import</translation>
    </message>
    <message>
        <source>[Remove] - remove the currently selected file(s) from the list</source>
        <translation>[Odobrať] - odoberie zo zoznamu vybraný/é súbory</translation>
    </message>
    <message>
        <source>[Remove All] - remove all the files in the list</source>
        <translation>[Odobrať všetko] - zo zoznamu odoberie všetky súbory</translation>
    </message>
    <message>
        <source>[SRID] - Reference ID for the shapefiles to be imported</source>
        <translation>[SRID] - Referenčné ID pre importované súbory shape</translation>
    </message>
    <message>
        <source>[Use Default (SRID)] - set SRID to -1</source>
        <translation>[Použiť predvolené (SRID)] - nastaví SRID na -1</translation>
    </message>
    <message>
        <source>[Geometry Column Name] - name of the geometry column in the database</source>
        <translation>[Meno stĺpca s geometriou] - meno stĺpca v databáze obsahujúceho geometriu</translation>
    </message>
    <message>
        <source>[Use Default (Geometry Column Name)] - set column name to &apos;the_geom&apos;</source>
        <translation>[Použiť predvolené (Meno stĺpca s geometriou)] - nastaví meno stĺpca na &apos;the_geom&apos;</translation>
    </message>
    <message>
        <source>[Glogal Schema] - set the schema for all files to be imported into</source>
        <translation>[Globálna schéma] - nastaví schému, do ktorej budú importované všetky súbory</translation>
    </message>
    <message>
        <source>[Import] - import the current shapefiles in the list</source>
        <translation>[Import] - naimportuje súbory shape z aktuálneho zoznamu</translation>
    </message>
    <message>
        <source>[Quit] - quit the program
</source>
        <translation>[Koniec] - ukončí program

</translation>
    </message>
    <message>
        <source>[Help] - display this help dialog</source>
        <translation>[Pomocník] - zobrazí tento pomocný text</translation>
    </message>
    <message>
        <source>Import Shapefiles</source>
        <translation>Importovať súbory shape</translation>
    </message>
    <message>
        <source>You need to specify a Connection first</source>
        <translation>Najprv je potrebné vybrať Spojenie</translation>
    </message>
    <message>
        <source>Connection failed - Check settings and try again</source>
        <translation>Spojenie zlyhalo - skontrolujte nastavenia a skúste znova</translation>
    </message>
    <message>
        <source>You need to add shapefiles to the list first</source>
        <translation>Najskôr je potrebné pridať do zoznamu súbory shape</translation>
    </message>
    <message>
        <source>Importing files</source>
        <translation>Importujú sa súbory</translation>
    </message>
    <message>
        <source>Cancel</source>
        <translation>Zrušiť</translation>
    </message>
    <message>
        <source>Progress</source>
        <translation>Priebeh</translation>
    </message>
    <message>
        <source>Problem inserting features from file:</source>
        <translation>Pri vkladaní objektov zo súboru nastal problém:</translation>
    </message>
    <message>
        <source>Invalid table name.</source>
        <translation>Nesprávny názov tabuľky.</translation>
    </message>
    <message>
        <source>No fields detected.</source>
        <translation>Nenašli sa žiadne polia.</translation>
    </message>
    <message>
        <source>The following fields are duplicates:</source>
        <translation>Nasledujúce polia sú zdvojené:</translation>
    </message>
    <message>
        <source>Import Shapefiles - Relation Exists</source>
        <translation>Import súborov shape - Relácia existuje</translation>
    </message>
    <message>
        <source>The Shapefile:</source>
        <translation>Súbor shape:</translation>
    </message>
    <message>
        <source>will use [</source>
        <translation>použije reláciu [</translation>
    </message>
    <message>
        <source>] relation for its data,</source>
        <translation>] pre svoje údaje,</translation>
    </message>
    <message>
        <source>which already exists and possibly contains data.</source>
        <translation>ktoré už existujú a pravdepodobne obsahujú údaje.</translation>
    </message>
    <message>
        <source>To avoid data loss change the &quot;DB Relation Name&quot;</source>
        <translation>Aby bolo možné vyhnúť sa strate údajov, je potrebné zmeniť &quot;DB relation Name&quot;</translation>
    </message>
    <message>
        <source>for this Shapefile in the main dialog file list.</source>
        <translation>pre tento súbor shape v hlavnom dialógu zoznamu súborov.</translation>
    </message>
    <message>
        <source>Do you want to overwrite the [</source>
        <translation>Chcete prepísať reláciu [</translation>
    </message>
    <message>
        <source>] relation?</source>
        <translation>]?</translation>
    </message>
    <message>
        <source>File Name</source>
        <translation>Meno súboru</translation>
    </message>
    <message>
        <source>Feature Class</source>
        <translation>Trieda objektu</translation>
    </message>
    <message>
        <source>Features</source>
        <translation>Objekty</translation>
    </message>
    <message>
        <source>DB Relation Name</source>
        <translation>Názov relácie v databáze</translation>
    </message>
    <message>
        <source>Schema</source>
        <translation>Schéma</translation>
    </message>
    <message>
        <source>Add Shapefiles</source>
        <translation>Pridať súbory Shape</translation>
    </message>
    <message>
        <source>Shapefiles (*.shp);;All files (*.*)</source>
        <translation>Súbory Shape (*.shp);;Všetky súbory (*.*)</translation>
    </message>
    <message>
        <source>PostGIS not available</source>
        <translation>PostGIS nie je dostupný</translation>
    </message>
    <message>
        <source>&lt;p&gt;The chosen database does not have PostGIS installed, but this is required for storage of spatial data.&lt;/p&gt;</source>
        <translation>&lt;p&gt;Vybraná databáza nemá nainštalované rozšírenie PostGIS, ktoré je nevyhnutné na ukladanie priestorových údajov.&lt;/p&gt;</translation>
    </message>
    <message>
        <source>Checking to see if </source>
        <translation type="obsolete">Zisťuje sa či </translation>
    </message>
    <message>
        <source>&lt;p&gt;Error while executing the SQL:&lt;/p&gt;&lt;p&gt;</source>
        <translation>&lt;p&gt;Chyba pri vykonávaní SQL dopytu:&lt;/p&gt;&lt;p&gt;</translation>
    </message>
    <message>
        <source>&lt;/p&gt;&lt;p&gt;The database said:</source>
        <translation>&lt;/p&gt;&lt;p&gt;Hlásenie z databázy:</translation>
    </message>
    <message>
        <source>%1 of %2 shapefiles could not be imported.</source>
        <translation>%1 z %2 súborov Shape sa nepodarilo naimportovať.</translation>
    </message>
    <message>
        <source>Password for </source>
        <translation>Heslo pre užívateľa</translation>
    </message>
    <message>
        <source>Please enter your password:</source>
        <translation>Prosím, zadajte vaše heslo:</translation>
    </message>
</context>
<context>
    <name>QgsSpitBase</name>
    <message>
        <source>SPIT - Shapefile to PostGIS Import Tool</source>
        <translation>SPIT - Nástroj na import súboru shape do PostGIS</translation>
    </message>
    <message>
        <source>PostgreSQL Connections</source>
        <translation>Spojenia PostgreSQL</translation>
    </message>
    <message>
        <source>Remove</source>
        <translation>Odobrať</translation>
    </message>
    <message>
        <source>Remove All</source>
        <translation>Odobrať všetko</translation>
    </message>
    <message>
        <source>Global Schema</source>
        <translation>Globálna schéma</translation>
    </message>
    <message>
        <source>Add</source>
        <translation>Pridať</translation>
    </message>
    <message>
        <source>Add a shapefile to the list of files to be imported</source>
        <translation>Pridať súbor shape do zoznamu importovaných súborov</translation>
    </message>
    <message>
        <source>Remove the selected shapefile from the import list</source>
        <translation>Odobrať označený súbor shape zo zoznamu importovaných súborov</translation>
    </message>
    <message>
        <source>Remove all the shapefiles from the import list</source>
        <translation>Odobrať všetky shape súbory zo zoznamu súborov určených na import</translation>
    </message>
    <message>
        <source>Set the SRID to the default value</source>
        <translation>Nastaví SRID na predvolenú hodnotu</translation>
    </message>
    <message>
        <source>Set the geometry column name to the default value</source>
        <translation>Nastaví názov stĺpca s geomteriou na predvolenú hodnotu</translation>
    </message>
    <message>
        <source>New</source>
        <translation>Nové</translation>
    </message>
    <message>
        <source>Create a new PostGIS connection</source>
        <translation>Vytvoriť spojenie PostGIS</translation>
    </message>
    <message>
        <source>Remove the current PostGIS connection</source>
        <translation>Odobrať aktuálne spojenie PostGIS</translation>
    </message>
    <message>
        <source>Connect</source>
        <translation>Pripojiť</translation>
    </message>
    <message>
        <source>Edit</source>
        <translation>Upraviť</translation>
    </message>
    <message>
        <source>Edit the current PostGIS connection</source>
        <translation>Upraviť toto spojenie PostGIS</translation>
    </message>
    <message>
        <source>Import options and shapefile list</source>
        <translation>Možnosti importu a zoznam súborov shape</translation>
    </message>
    <message>
        <source>Use Default SRID or specify here</source>
        <translation>Použiť predvolené SRID alebo ho zadať tu</translation>
    </message>
    <message>
        <source>Use Default Geometry Column Name or specify here</source>
        <translation>Použiť predvolené meno stĺpca s geometriou alebo ho zadať tu</translation>
    </message>
    <message>
        <source>Primary Key Column Name</source>
        <translation>Meno stĺpca s primárnym kľúčom</translation>
    </message>
    <message>
        <source>Connect to PostGIS</source>
        <translation>Pripojí sa k PostGIS</translation>
    </message>
</context>
<context>
    <name>QgsSpitPlugin</name>
    <message>
        <source>&amp;Import Shapefiles to PostgreSQL</source>
        <translation>&amp;Import súborov Shape do PostgreSQL</translation>
    </message>
    <message>
        <source>Import shapefiles into a PostGIS-enabled PostgreSQL database. The schema and field names can be customized on import</source>
        <translation>Importuje shape súbory do databázy PostgreSQL s rozšírením PostGIS. Schému a názvy polí možno pred importom upraviť</translation>
    </message>
    <message>
        <source>&amp;Spit</source>
        <translation>&amp;Spit</translation>
    </message>
</context>
<context>
    <name>QgsTINInterpolatorDialog</name>
    <message>
        <source>Linear interpolation</source>
        <translation>Lineárna interpolácia</translation>
    </message>
</context>
<context>
    <name>QgsTINInterpolatorDialogBase</name>
    <message>
        <source>Triangle based interpolation</source>
        <translation>Interpolácia založená na trojuholníkoch</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:12pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;This interpolator provides different methods for interpolation in a triangular irregular network (TIN).&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation type="unfinished">&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;Sans Serif&apos;; font-size:12pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;&quot;&gt;Tento modul zabezpečuje rôzne metódy pre interpoláciu v nepravidelnej sieti trojuholníkov (TIN).&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>Interpolation method:</source>
        <translation>Interpolačná metóda:</translation>
    </message>
</context>
<context>
    <name>QgsUniqueValueDialog</name>
    <message>
        <source>Confirm Delete</source>
        <translation>Potvrdenie mazania</translation>
    </message>
    <message>
        <source>The classification field was changed from &apos;%1&apos; to &apos;%2&apos;.
Should the existing classes be deleted before classification?</source>
        <translation>Klasifikačné pole sa zmenilo z %1 na %2.
Majú byť pred klasifikáciou zmazané dosiaľ existujúce triedy?</translation>
    </message>
</context>
<context>
    <name>QgsUniqueValueDialogBase</name>
    <message>
        <source>Form1</source>
        <translation>Jedinečná hodnota</translation>
    </message>
    <message>
        <source>Classify</source>
        <translation>Klasifikovať</translation>
    </message>
    <message>
        <source>Classification field</source>
        <translation>Klasifikačné pole</translation>
    </message>
    <message>
        <source>Add class</source>
        <translation>Pridať triedu</translation>
    </message>
    <message>
        <source>Delete classes</source>
        <translation>Zmazať triedy</translation>
    </message>
    <message>
        <source>Randomize Colors</source>
        <translation>Náhodné farby</translation>
    </message>
    <message>
        <source>Reset Colors</source>
        <translation>Znovunastaviť farby</translation>
    </message>
</context>
<context>
    <name>QgsVectorLayer</name>
    <message>
        <source>Could not commit the added features.</source>
        <translation type="obsolete">Nemožno zapísať pridané objekty.</translation>
    </message>
    <message>
        <source>No other types of changes will be committed at this time.</source>
        <translation type="obsolete">Žiadne iné typy zmien nebudú zapísané.</translation>
    </message>
    <message>
        <source>Could not commit the changed attributes.</source>
        <translation type="obsolete">Nemožno zapísať zmeny atribútov.</translation>
    </message>
    <message>
        <source>However, the added features were committed OK.</source>
        <translation type="obsolete">Avšak pridané objekty boli zapísané v poriadku.</translation>
    </message>
    <message>
        <source>Could not commit the changed geometries.</source>
        <translation type="obsolete">Nemožno zapísať zmenené geometrie.</translation>
    </message>
    <message>
        <source>However, the changed attributes were committed OK.</source>
        <translation type="obsolete">Avšak zmenené atribúty boli zapísané v poriadku.</translation>
    </message>
    <message>
        <source>Could not commit the deleted features.</source>
        <translation type="obsolete">Nemožno zapísať vymazanie objektov.</translation>
    </message>
    <message>
        <source>However, the changed geometries were committed OK.</source>
        <translation type="obsolete">Avšak zmeny geometrií boli zapísané v poriadku.</translation>
    </message>
    <message>
        <source>ERROR: no provider</source>
        <translation>CHYBA: žiadny správca</translation>
    </message>
    <message>
        <source>ERROR: layer not editable</source>
        <translation>CHYBA: vrstva nie je upravovateľná</translation>
    </message>
    <message>
        <source>SUCCESS: %1 attributes added.</source>
        <translation type="unfinished">ÚSPECH: %1 atribútov pridaných.</translation>
    </message>
    <message>
        <source>ERROR: %1 new attributes not added</source>
        <translation type="unfinished">CHYBA: %1 nových atribútov nebolo pridaných</translation>
    </message>
    <message>
        <source>SUCCESS: %1 attributes deleted.</source>
        <translation type="unfinished">ÚSPECH: %1 atribútov vymazaných.</translation>
    </message>
    <message>
        <source>ERROR: %1 attributes not deleted.</source>
        <translation type="unfinished">CHYBA: %1 atribútov nebolo vymazaných.</translation>
    </message>
    <message>
        <source>SUCCESS: attribute %1 was added.</source>
        <translation type="unfinished">ÚSPECH: atribút %1 bol pridaný.</translation>
    </message>
    <message>
        <source>ERROR: attribute %1 not added</source>
        <translation type="unfinished">CHYBA: atribút %1 nebol pridaný</translation>
    </message>
    <message>
        <source>SUCCESS: %1 attribute values changed.</source>
        <translation type="unfinished">ÚSPECH: %1 atribútových hodnôt bolo zmenených.</translation>
    </message>
    <message>
        <source>ERROR: %1 attribute value changes not applied.</source>
        <translation type="unfinished">CHYBA: %1 zmien atribútových hodnôt nebolo uskutočnených.</translation>
    </message>
    <message>
        <source>SUCCESS: %1 features added.</source>
        <translation type="unfinished">ÚSPECH: %1 objektov pridaných.</translation>
    </message>
    <message>
        <source>ERROR: %1 features not added.</source>
        <translation type="unfinished">CHYBA: %1 objektov nebolo pridaných.</translation>
    </message>
    <message>
        <source>SUCCESS: %1 geometries were changed.</source>
        <translation type="unfinished">ÚSPECH: %1 geometrií bolo zmenených.</translation>
    </message>
    <message>
        <source>ERROR: %1 geometries not changed.</source>
        <translation type="unfinished">CHYBA: %1 geometrií nezmenených.</translation>
    </message>
    <message>
        <source>SUCCESS: %1 features deleted.</source>
        <translation type="unfinished">ÚSPECH: %1 objektov vymazaných.</translation>
    </message>
    <message>
        <source>ERROR: %1 features not deleted.</source>
        <translation type="unfinished">CHYBA: %1 objektov nebolo vymazaných.</translation>
    </message>
    <message>
        <source>No renderer object</source>
        <translation type="unfinished">Žiadny objekt vykresľovača</translation>
    </message>
    <message>
        <source>Classification field not found</source>
        <translation type="unfinished">Nenašlo sa klasifikačné pole</translation>
    </message>
</context>
<context>
    <name>QgsVectorLayerProperties</name>
    <message>
        <source>Transparency: </source>
        <translation>Priehľadnosť: </translation>
    </message>
    <message>
        <source>Single Symbol</source>
        <translation>Jeden symbol</translation>
    </message>
    <message>
        <source>Graduated Symbol</source>
        <translation>Odstupňovaný symbol</translation>
    </message>
    <message>
        <source>Continuous Color</source>
        <translation>Spojitá farba</translation>
    </message>
    <message>
        <source>Unique Value</source>
        <translation>Jedinečná hodnota</translation>
    </message>
    <message>
        <source>This button opens the PostgreSQL query builder and allows you to create a subset of features to display on the map canvas rather than displaying all features in the layer</source>
        <translation type="unfinished">Týmto tlačidlom sa otvára nástroj na tvorbu dopytov pre PostgreSQL, čo umožňuje vytvoriť na rozdiel od zobrazenia všetkých objektov vrstvy, len určitú podskupinu objektov, ktoré sa zobrazia na mapovom plátne </translation>
    </message>
    <message>
        <source>The query used to limit the features in the layer is shown here. This is currently only supported for PostgreSQL layers. To enter or modify the query, click on the Query Builder button</source>
        <translation type="unfinished">Tu sa zobrazuje dopyt použitý na vymedzenie objektov vo vrstve. V súčasnosti sú dopyty podporované len pre vrstvy PostgreSQL. Zadávať, alebo upravovať dopyt možno kliknutím na tlačidlo Tvorba dopytov</translation>
    </message>
    <message>
        <source>Spatial Index</source>
        <translation>Priestorový index</translation>
    </message>
    <message>
        <source>Creation of spatial index successfull</source>
        <translation type="obsolete">Vytvorenie priestorového indexu bolo úspešné</translation>
    </message>
    <message>
        <source>Creation of spatial index failed</source>
        <translation>Vytvorenie priestorového indexu zlyhalo</translation>
    </message>
    <message>
        <source>General:</source>
        <translation>Všeobecné:</translation>
    </message>
    <message>
        <source>Storage type of this layer : </source>
        <translation> Typ uloženia tejto vrstvy:</translation>
    </message>
    <message>
        <source>Geometry type of the features in this layer : </source>
        <translation>Geometrický typ objektov v tejto vrstve: </translation>
    </message>
    <message>
        <source>The number of features in this layer : </source>
        <translation>Počet objektov v tejto vrstve: </translation>
    </message>
    <message>
        <source>Editing capabilities of this layer : </source>
        <translation> Možnosti úprav tejto vrstvy: </translation>
    </message>
    <message>
        <source>Extents:</source>
        <translation>Rozsah:</translation>
    </message>
    <message>
        <source>In layer spatial reference system units : </source>
        <translation> V jednotkách referenčného priestorového systému vrstvy: </translation>
    </message>
    <message>
        <source>xMin,yMin </source>
        <translation>xMin,yMin </translation>
    </message>
    <message>
        <source> : xMax,yMax </source>
        <translation> : xMax,yMax </translation>
    </message>
    <message>
        <source>In project spatial reference system units : </source>
        <translation>V jednotkách priestorového referenčného systému projektu: </translation>
    </message>
    <message>
        <source>Layer Spatial Reference System:</source>
        <translation>Priestorový referenčný systém vrstvy:</translation>
    </message>
    <message>
        <source>Attribute field info:</source>
        <translation>Informácia o atribútovom poli:</translation>
    </message>
    <message>
        <source>Field</source>
        <translation>Pole</translation>
    </message>
    <message>
        <source>Type</source>
        <translation>Typ</translation>
    </message>
    <message>
        <source>Length</source>
        <translation>Dĺžka</translation>
    </message>
    <message>
        <source>Precision</source>
        <translation>Presnosť</translation>
    </message>
    <message>
        <source>Source for this layer : </source>
        <translation> Zdroj pre túto vrstvu:</translation>
    </message>
    <message>
        <source>Layer comment: </source>
        <translation>Komentáre k vrstve: </translation>
    </message>
    <message>
        <source>Comment</source>
        <translation>Komentár</translation>
    </message>
    <message>
        <source>Default Style</source>
        <translation>Predvolený štýl</translation>
    </message>
    <message>
        <source>QGIS Layer Style File (*.qml)</source>
        <translation>Súbor so štýlom vrstvy pre QGIS (.qml)</translation>
    </message>
    <message>
        <source>QGIS</source>
        <translation>QGIS</translation>
    </message>
    <message>
        <source>Unknown style format: </source>
        <translation>Neznámy formát stýlu: </translation>
    </message>
    <message>
        <source>Saved Style</source>
        <translation>Uložený štýl</translation>
    </message>
    <message>
        <source>id</source>
        <translation>id</translation>
    </message>
    <message>
        <source>name</source>
        <translation>meno</translation>
    </message>
    <message>
        <source>type</source>
        <translation>typ</translation>
    </message>
    <message>
        <source>length</source>
        <translation>dĺžka</translation>
    </message>
    <message>
        <source>precision</source>
        <translation>presnosť</translation>
    </message>
    <message>
        <source>comment</source>
        <translation>komentár</translation>
    </message>
    <message>
        <source>edit widget</source>
        <translation type="unfinished">upraviť widget</translation>
    </message>
    <message>
        <source>values</source>
        <translation>hodnoty</translation>
    </message>
    <message>
        <source>line edit</source>
        <translation type="unfinished">úprava línie</translation>
    </message>
    <message>
        <source>unique values</source>
        <translation type="unfinished">jedinečné hodnoty</translation>
    </message>
    <message>
        <source>unique values (editable)</source>
        <translation type="unfinished">jedinečné hodnoty (upravovateľné)</translation>
    </message>
    <message>
        <source>value map</source>
        <translation type="unfinished">mapa hodnôt</translation>
    </message>
    <message>
        <source>classification</source>
        <translation type="unfinished">klasifikácia</translation>
    </message>
    <message>
        <source>range (editable)</source>
        <translation type="unfinished">rozsah (upravovateľný)</translation>
    </message>
    <message>
        <source>range (slider)</source>
        <translation type="unfinished">rozsah (posuvník)</translation>
    </message>
    <message>
        <source>file name</source>
        <translation>meno súboru</translation>
    </message>
    <message>
        <source>Name conflict</source>
        <translation>Konflikt názov</translation>
    </message>
    <message>
        <source>The attribute could not be inserted. The name already exists in the table.</source>
        <translation>Tento atribút nemožno vložiť. Atribút s takýmto názvom už v tabuľke existuje.</translation>
    </message>
    <message>
        <source>Creation of spatial index successful</source>
        <translation>Vytvorenie priestorového indexu bolo úspešné</translation>
    </message>
</context>
<context>
    <name>QgsVectorLayerPropertiesBase</name>
    <message>
        <source>Layer Properties</source>
        <translation>Vlastnosti vrstvy</translation>
    </message>
    <message>
        <source>Symbology</source>
        <translation>Symbolika</translation>
    </message>
    <message>
        <source>General</source>
        <translation>Všeobecné</translation>
    </message>
    <message>
        <source>Use scale dependent rendering</source>
        <translation>Používať vykresľovanie v závislosti od mierky</translation>
    </message>
    <message>
        <source>Minimum scale at which this layer will be displayed. </source>
        <translation>Minimálna mierka pri ktorej bude táto vrstva zobrazená.</translation>
    </message>
    <message>
        <source>Maximum scale at which this layer will be displayed. </source>
        <translation>Maximálna mierka pri ktorej bude táto vrstva zobrazená.</translation>
    </message>
    <message>
        <source>Display name</source>
        <translation>Zobraziť meno</translation>
    </message>
    <message>
        <source>Use this control to set which field is placed at the top level of the Identify Results dialog box.</source>
        <translation>Používa sa na nastavenie poľa, ktoré bude umiestnené na najvyššej úrovni dialógového okna Výsledky identifikácie.</translation>
    </message>
    <message>
        <source>Display field for the Identify Results dialog box</source>
        <translation>Pole zobrazované v dialógovom okne Výsledky identifikácie</translation>
    </message>
    <message>
        <source>This sets the display field for the Identify Results dialog box</source>
        <translation>Týmto sa nastaví pole zobrazované v dialógovom okne Výsledky identifikácie</translation>
    </message>
    <message>
        <source>Display field</source>
        <translation>Zobrazované pole</translation>
    </message>
    <message>
        <source>Subset</source>
        <translation>Podskupina</translation>
    </message>
    <message>
        <source>Query Builder</source>
        <translation>Tvorba dopytov</translation>
    </message>
    <message>
        <source>Create Spatial Index</source>
        <translation>Vytvoriť priestorový index</translation>
    </message>
    <message>
        <source>Metadata</source>
        <translation>Meta údaje</translation>
    </message>
    <message>
        <source>Labels</source>
        <translation>Popisy</translation>
    </message>
    <message>
        <source>Display labels</source>
        <translation>Zobraziť popisy</translation>
    </message>
    <message>
        <source>Actions</source>
        <translation>Akcie</translation>
    </message>
    <message>
        <source>Restore Default Style</source>
        <translation>Obnoviť predvolený štýl</translation>
    </message>
    <message>
        <source>Save As Default</source>
        <translation>Uložiť ako predvolený</translation>
    </message>
    <message>
        <source>Load Style ...</source>
        <translation>Nahrať štýl ...</translation>
    </message>
    <message>
        <source>Save Style ...</source>
        <translation>Uložiť štýl ...</translation>
    </message>
    <message>
        <source>Legend type</source>
        <translation>Druh legendy</translation>
    </message>
    <message>
        <source>Transparency</source>
        <translation>Priehľadnosť</translation>
    </message>
    <message>
        <source>Options</source>
        <translation>Nastavenia</translation>
    </message>
    <message>
        <source>Change SRS</source>
        <translation type="obsolete">Zmeniť SRS</translation>
    </message>
    <message>
        <source>Maximum</source>
        <translation>Maximum</translation>
    </message>
    <message>
        <source>Minimum</source>
        <translation>Minimum</translation>
    </message>
    <message>
        <source>Change CRS</source>
        <translation type="unfinished">Zmeniť CRS</translation>
    </message>
    <message>
        <source>Attributes</source>
        <translation>Atribúty</translation>
    </message>
    <message>
        <source>New column</source>
        <translation>Nový stĺpec</translation>
    </message>
    <message>
        <source>Ctrl+N</source>
        <translation>Ctrl+N</translation>
    </message>
    <message>
        <source>Delete column</source>
        <translation>Vymazať stĺpec</translation>
    </message>
    <message>
        <source>Ctrl+X</source>
        <translation>Ctrl+X</translation>
    </message>
    <message>
        <source>Toggle editing mode</source>
        <translation>Prepnúť do režimu úprav</translation>
    </message>
    <message>
        <source>Click to toggle table editing</source>
        <translation>Po kliknutí sa prepne tabuľka do režimu úprav</translation>
    </message>
</context>
<context>
    <name>QgsVectorSymbologyWidgetBase</name>
    <message>
        <source>Form2</source>
        <translation type="obsolete">Symbolika vektorovej vrstvy</translation>
    </message>
    <message>
        <source>Label</source>
        <translation type="obsolete">Popis</translation>
    </message>
    <message>
        <source>Min</source>
        <translation type="obsolete">Min</translation>
    </message>
    <message>
        <source>Max</source>
        <translation type="obsolete">Max</translation>
    </message>
    <message>
        <source>Symbol Classes:</source>
        <translation type="obsolete">Triedy symbolov:</translation>
    </message>
    <message>
        <source>Count:</source>
        <translation type="obsolete">Počet:</translation>
    </message>
    <message>
        <source>Mode:</source>
        <translation type="obsolete">Mód:</translation>
    </message>
    <message>
        <source>Field:</source>
        <translation type="obsolete">Pole:</translation>
    </message>
</context>
<context>
    <name>QgsWFSPlugin</name>
    <message>
        <source>&amp;Add WFS layer</source>
        <translation>&amp;Pridať vrstvu WFS</translation>
    </message>
</context>
<context>
    <name>QgsWFSProvider</name>
    <message>
        <source>unknown</source>
        <translation type="unfinished"> neznáma veľkosť</translation>
    </message>
    <message>
        <source>received %1 bytes from %2</source>
        <translation type="unfinished">prijatých %1 z %2 bajtov</translation>
    </message>
</context>
<context>
    <name>QgsWFSSourceSelect</name>
    <message>
        <source>Are you sure you want to remove the </source>
        <translation>Ste si istý, že chcete odstrániť spojenie </translation>
    </message>
    <message>
        <source> connection and all associated settings?</source>
        <translation>  a všetky s ním súvisiace nastavenia?</translation>
    </message>
    <message>
        <source>Confirm Delete</source>
        <translation>Potvrdenie mazania</translation>
    </message>
</context>
<context>
    <name>QgsWFSSourceSelectBase</name>
    <message>
        <source>Server Connections</source>
        <translation>Spojenia</translation>
    </message>
    <message>
        <source>&amp;New</source>
        <translation>&amp;Nové</translation>
    </message>
    <message>
        <source>Delete</source>
        <translation>Vymazať</translation>
    </message>
    <message>
        <source>Edit</source>
        <translation>Upraviť</translation>
    </message>
    <message>
        <source>C&amp;onnect</source>
        <translation>&amp;Pripojiť</translation>
    </message>
    <message>
        <source>Coordinate Reference System</source>
        <translation>Referenčný súradnicový systém</translation>
    </message>
    <message>
        <source>Change ...</source>
        <translation>Zmeniť ...</translation>
    </message>
    <message>
        <source>Title</source>
        <translation>Titulok</translation>
    </message>
    <message>
        <source>Name</source>
        <translation>Názov</translation>
    </message>
    <message>
        <source>Abstract</source>
        <translation>Abstrakt</translation>
    </message>
    <message>
        <source>Add WFS Layer from a Server</source>
        <translation>Pridanie vrstvy WFS (zo servera)</translation>
    </message>
</context>
<context>
    <name>QgsWmsProvider</name>
    <message>
        <source>Server Properties:</source>
        <translation>Vlastnosti servera:</translation>
    </message>
    <message>
        <source>Property</source>
        <translation>Vlastnosť</translation>
    </message>
    <message>
        <source>Value</source>
        <translation>Hodnota</translation>
    </message>
    <message>
        <source>WMS Version</source>
        <translation>Verzia WMS</translation>
    </message>
    <message>
        <source>Title</source>
        <translation>Názov</translation>
    </message>
    <message>
        <source>Abstract</source>
        <translation>Abstrakt</translation>
    </message>
    <message>
        <source>Keywords</source>
        <translation>Kľúčové slová</translation>
    </message>
    <message>
        <source>Online Resource</source>
        <translation>Online zdroje</translation>
    </message>
    <message>
        <source>Contact Person</source>
        <translation>Kontaktná osoba</translation>
    </message>
    <message>
        <source>Fees</source>
        <translation>Poplatky</translation>
    </message>
    <message>
        <source>Access Constraints</source>
        <translation>Obmedzenia prístupu</translation>
    </message>
    <message>
        <source>Image Formats</source>
        <translation>Formát obrázkov</translation>
    </message>
    <message>
        <source>Layer Count</source>
        <translation>Počet vrstiev</translation>
    </message>
    <message>
        <source>Layer Properties: </source>
        <translation>Vlastnosti vrstvy: </translation>
    </message>
    <message>
        <source>Selected</source>
        <translation>Vybrané</translation>
    </message>
    <message>
        <source>Yes</source>
        <translation>Áno</translation>
    </message>
    <message>
        <source>No</source>
        <translation>Nie</translation>
    </message>
    <message>
        <source>Visibility</source>
        <translation>Viditeľnosť</translation>
    </message>
    <message>
        <source>Visible</source>
        <translation>Viditeľná</translation>
    </message>
    <message>
        <source>Hidden</source>
        <translation>Skrytá</translation>
    </message>
    <message>
        <source>n/a</source>
        <translation>--</translation>
    </message>
    <message>
        <source>Available in CRS</source>
        <translation type="obsolete">Dostupná v CRS</translation>
    </message>
    <message>
        <source>Available in style</source>
        <translation>Dostupná v štýle</translation>
    </message>
    <message>
        <source>Name</source>
        <translation>Meno</translation>
    </message>
    <message>
        <source>HTTP Exception</source>
        <translation>HTTP výnimka</translation>
    </message>
    <message>
        <source>Tried URL: </source>
        <translation>Použité URL: </translation>
    </message>
    <message>
        <source>WMS Service Exception</source>
        <translation>Výnimka služby WMS</translation>
    </message>
    <message>
        <source>DOM Exception</source>
        <translation type="obsolete">DOM výnimka</translation>
    </message>
    <message>
        <source>Could not get WMS capabilities: %1 at line %2 column %3</source>
        <translation>Nemožno získať WMS službu: %1 na riadku %2 stĺpec %3</translation>
    </message>
    <message>
        <source>This is probably due to an incorrect WMS Server URL.</source>
        <translation>Je to spôsobené pravdepodobne nesprávnym URL servera WMS.</translation>
    </message>
    <message>
        <source>Could not get WMS capabilities in the expected format (DTD): no %1 or %2 found</source>
        <translation>Nemožno získať WMS službu v očakávanom formáte (DTD): žiadny %1 alebo %2</translation>
    </message>
    <message>
        <source>Could not get WMS Service Exception at %1: %2 at line %3 column %4</source>
        <translation>Nemožno získať WMS službu Výnimka na %1: %2 v riadku %3 stĺpec %4</translation>
    </message>
    <message>
        <source>Request contains a Format not offered by the server.</source>
        <translation>Požiadavka obsahuje Format, ktorý tento server neponúka.</translation>
    </message>
    <message>
        <source>Request contains a CRS not offered by the server for one or more of the Layers in the request.</source>
        <translation>Požiadavka obsahuje CRS (referenčný súradnicový systém) ktorý nie je v ponuke servera pre jednu alebo viac vrstiev v tejto požiadavke.</translation>
    </message>
    <message>
        <source>Request contains a SRS not offered by the server for one or more of the Layers in the request.</source>
        <translation>Požiadavka obsahuje priestorový referenčný systém, ktorý nie je v ponuke servera pre jednu alebo viac vrstiev v požiadavke.</translation>
    </message>
    <message>
        <source>GetMap request is for a Layer not offered by the server, or GetFeatureInfo request is for a Layer not shown on the map.</source>
        <translation>Požiadavka GetMap je pre vrstvu, ktorú server neponúka, alebo požiadavka GetFeatureInfo je pre vrstvu, ktorá nie je na mape.</translation>
    </message>
    <message>
        <source>Request is for a Layer in a Style not offered by the server.</source>
        <translation>Požiadavka je na vrstvu v štýle ktorý neponúka tento server.</translation>
    </message>
    <message>
        <source>GetFeatureInfo request is applied to a Layer which is not declared queryable.</source>
        <translation>Požiadavka GetFeatureInfo je použitá na vrstvu, ktorá nie je uvedená ako dopytovateľná.</translation>
    </message>
    <message>
        <source>GetFeatureInfo request contains invalid X or Y value.</source>
        <translation>Požiadavka GetFeatureInfo obsahuje chybné hodnoty X alebo Y.</translation>
    </message>
    <message>
        <source>Value of (optional) UpdateSequence parameter in GetCapabilities request is equal to current value of service metadata update sequence number.</source>
        <translation>Hodnota (voliteľného) parametra UpdateSequence v požiadavke GetCapabilities je rovnaká s aktuálnou hodnotou service metadata update sequence number.</translation>
    </message>
    <message>
        <source>Value of (optional) UpdateSequence parameter in GetCapabilities request is greater than current value of service metadata update sequence number.</source>
        <translation>Hodnota (voliteľného) parametra UpdateSequence v požiadavke GetCapabilities je väčšia ako aktuálna hodnota service metadata update sequence number.</translation>
    </message>
    <message>
        <source>Request does not include a sample dimension value, and the server did not declare a default value for that dimension.</source>
        <translation>Požiadavka neobsahuje vzorku hodnoty veľkosti a server neurčil predvolenú hodnotu pre takýto rozmer.</translation>
    </message>
    <message>
        <source>Request contains an invalid sample dimension value.</source>
        <translation>Požiadavka obsahuje neplatnú vzorku hodnoty rozmeru.</translation>
    </message>
    <message>
        <source>Request is for an optional operation that is not supported by the server.</source>
        <translation>Požiadavka smeruje na voliteľnú operáciu, ktorú tento server nepodporuje.</translation>
    </message>
    <message>
        <source>(Unknown error code from a post-1.3 WMS server)</source>
        <translation>(Neznámy chybový kód z post-1.3 WMS servera)</translation>
    </message>
    <message>
        <source>The WMS vendor also reported: </source>
        <translation> Správa od poskytovateľa WMS: </translation>
    </message>
    <message>
        <source>This is probably due to a bug in the QGIS program.  Please report this error.</source>
        <translation type="obsolete">  Je to pravdepodobne spôsobené chybou v programe QGIS. Prosím oznámte túto chybu.</translation>
    </message>
    <message>
        <source>Identify Formats</source>
        <translation type="unfinished">Identifikácia formátov</translation>
    </message>
    <message>
        <source>Can Identify</source>
        <translation>Možno identifikovať</translation>
    </message>
    <message>
        <source>Can be Transparent</source>
        <translation>Môže byť priehľadná</translation>
    </message>
    <message>
        <source>Can Zoom In</source>
        <translation>Možno priblížiť</translation>
    </message>
    <message>
        <source>Cascade Count</source>
        <translation type="unfinished">Kaskádové počítanie (počet)</translation>
    </message>
    <message>
        <source>Fixed Width</source>
        <translation>Pevná sírka</translation>
    </message>
    <message>
        <source>Fixed Height</source>
        <translation>Pevná výška</translation>
    </message>
    <message>
        <source>WGS 84 Bounding Box</source>
        <translation>WGS 84 vymedzujúci/ohraničujúci obdĺžnik</translation>
    </message>
    <message>
        <source>Layer cannot be queried.</source>
        <translation>Na vrstvu nemožno robiť dopyty.</translation>
    </message>
    <message>
        <source>Dom Exception</source>
        <translation>Výnimka DOM</translation>
    </message>
    <message>
        <source>Could not determine URL for GetMap from the WMS capabilities response</source>
        <translation type="obsolete">Z odpovede WMS capabilities nemožno zistiť URL pre GetMap</translation>
    </message>
</context>
<context>
    <name>QuickPrintGui</name>
    <message>
        <source>Portable Document Format (*.pdf)</source>
        <translation>Portable Document Format (*.pdf)</translation>
    </message>
    <message>
        <source>quickprint</source>
        <translation type="unfinished">Rýchlotlač</translation>
    </message>
    <message>
        <source>Unknown format: </source>
        <translation>Neznámy formát: </translation>
    </message>
</context>
<context>
    <name>QuickPrintGuiBase</name>
    <message>
        <source>QGIS Quick Print Plugin</source>
        <translation>Zásuvný modul QGISu na rýchlu tlač</translation>
    </message>
    <message>
        <source>Quick Print</source>
        <translation>Rýchla tlač</translation>
    </message>
    <message>
        <source>Map Title e.g. ACME inc.</source>
        <translation type="unfinished">Titulok mapy (napr. ACME inc.)</translation>
    </message>
    <message>
        <source>Map Name e.g. Water Features</source>
        <translation>Názov mapy (napr. Vodstvo)</translation>
    </message>
    <message>
        <source>Copyright</source>
        <translation>Autorské práva</translation>
    </message>
    <message>
        <source>Output</source>
        <translation>Výstup</translation>
    </message>
    <message>
        <source>Use last filename but incremented.</source>
        <translation>Použiť posledne použité meno súboru zväčšené o 1.</translation>
    </message>
    <message>
        <source>last used filename but incremented will be shown here</source>
        <translation>tu sa objaví naposledy použité meno súboru zväčšené o 1</translation>
    </message>
    <message>
        <source>Prompt for file name</source>
        <translation>Opýtať sa meno súboru</translation>
    </message>
    <message>
        <source>Note: If you want more control over the map layout please use the map composer function in QGIS.</source>
        <translation>Poznámka: Pokiaľ si želáte viac ovládacích prvkov nad mapovým výstupom použite v QGise funkciu Skladateľ máp.</translation>
    </message>
    <message>
        <source>Page Size</source>
        <translation>Veľkosť strany</translation>
    </message>
</context>
<context>
    <name>QuickPrintPlugin</name>
    <message>
        <source>Quick Print</source>
        <translation>Rýchla tlač</translation>
    </message>
    <message>
        <source>Replace this with a short description of the what the plugin does</source>
        <translation type="obsolete">Zameňte tento text s opisom činnosti daného zásuvného modulu</translation>
    </message>
    <message>
        <source>&amp;Quick Print</source>
        <translation>&amp;Rýchla tlač</translation>
    </message>
    <message>
        <source>Provides a qay to quickly produce a map with minimal user input.</source>
        <translation type="obsolete">Poskytuje spôsob ako rýchlo vytvoriť mapu s minimálnym vstupom od používateľa.</translation>
    </message>
    <message>
        <source>Provides a way to quickly produce a map with minimal user input.</source>
        <translation>Poskytuje spôsob ako rýchlo vytvoriť mapu s minimálnym vstupom od používateľa.</translation>
    </message>
</context>
<context>
    <name>RepositoryDetailsDialog</name>
    <message>
        <source>Repository details</source>
        <translation type="obsolete">Podrobnosti repozitára</translation>
    </message>
    <message>
        <source>Name:</source>
        <translation type="obsolete">Meno:</translation>
    </message>
    <message>
        <source>URL:</source>
        <translation type="obsolete">URL:</translation>
    </message>
    <message>
        <source>http://</source>
        <translation type="obsolete">http://</translation>
    </message>
</context>
<context>
    <name>[pluginname]GuiBase</name>
    <message>
        <source>QGIS Plugin Template</source>
        <translation>QGIS šablóna zásuvného modulu</translation>
    </message>
    <message>
        <source>Plugin Template</source>
        <translation>Šablóna zásuvného modulu</translation>
    </message>
</context>
<context>
    <name>dxf2shpConverter</name>
    <message>
        <source>Converts DXF files in Shapefile format</source>
        <translation>Prevedie súbory DXF do formátu Shape</translation>
    </message>
    <message>
        <source>&amp;Dxf2Shp</source>
        <translation>&amp;Dxf2Shp</translation>
    </message>
</context>
<context>
    <name>dxf2shpConverterGui</name>
    <message>
        <source>Fields description:
* Input DXF file: path to the DXF file to be converted
* Output Shp file: desired name of the shape file to be created
* Shp output file type: specifies the type of the output shape file
* Export text labels checkbox: if checked, an additional shp points layer will be created,   and the associated dbf table will contain informations about the &quot;TEXT&quot; fields found in the dxf file, and the text strings themselves

---
Developed by Paolo L. Scala, Barbara Rita Barricelli, Marco Padula
CNR, Milan Unit (Information Technology), Construction Technologies Institute.
For support send a mail to scala@itc.cnr.it
</source>
        <translation type="obsolete">Popis polí:
* Vstupný súbor DXF: cesta k súboru DXF určenému na prevod
* Výstupný súbor Shp: meno súboru ktorý bude vytvorený
* Typ výstupného súboru Shp: určuje typ výstupného súboru Shape
* Zaškrtávacie políčko &quot;Exportovať textové popisy&quot;: pokiaľ je zaškrtnuté, bude vytvorená ďalšia bodová vrstva vo formáte shp a k nej prislúchajúca tabuľka dbf bude obsahovať informácie o poliach &quot;TEXT&quot; nájdených v súbore dxf a im zodpovedajúce textové reťazce

---
Vyvinuli: Paolo L. Scala, Barbara Rita Barricelli, Marco Padula
CNR, Milan Unit (Information Technology), Construction Technologies Institute.
Podporu môžete získať zaslaním e-mailu na scala@itc.cnr.it</translation>
    </message>
    <message>
        <source>Choose a DXF file to open</source>
        <translation type="obsolete">Výber súboru DXF</translation>
    </message>
    <message>
        <source>Choose a file name to save to</source>
        <translation type="obsolete">Výber cieľového súboru</translation>
    </message>
    <message>
        <source>Dxf Importer</source>
        <translation>Importér súborov Dxf</translation>
    </message>
    <message>
        <source>Input Dxf file</source>
        <translation>Vstupný súbor Dxf</translation>
    </message>
    <message>
        <source>...</source>
        <translation>...</translation>
    </message>
    <message>
        <source>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;DejaVu Sans&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;MS Shell Dlg 2&apos;; font-size:8pt;&quot;&gt;&lt;span style=&quot; font-size:10pt;&quot;&gt;Output file&lt;/span&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</source>
        <translation>&lt;html&gt;&lt;head&gt;&lt;meta name=&quot;qrichtext&quot; content=&quot;1&quot; /&gt;&lt;style type=&quot;text/css&quot;&gt;
p, li { white-space: pre-wrap; }
&lt;/style&gt;&lt;/head&gt;&lt;body style=&quot; font-family:&apos;DejaVu Sans&apos;; font-size:10pt; font-weight:400; font-style:normal;&quot;&gt;
&lt;p style=&quot; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:&apos;MS Shell Dlg 2&apos;; font-size:8pt;&quot;&gt;&lt;span style=&quot; font-size:10pt;&quot;&gt;Výstupný súbor&lt;/span&gt;&lt;/p&gt;&lt;/body&gt;&lt;/html&gt;</translation>
    </message>
    <message>
        <source>Output file type</source>
        <translation>Typ výstupného súboru</translation>
    </message>
    <message>
        <source>Polyline</source>
        <translation>línia</translation>
    </message>
    <message>
        <source>Polygon</source>
        <translation>polygón</translation>
    </message>
    <message>
        <source>Point</source>
        <translation>bod</translation>
    </message>
    <message>
        <source>Export text labels</source>
        <translation>Exportovať textové popisy</translation>
    </message>
</context>
<context>
    <name>pluginname</name>
    <message>
        <source>Replace this with a short description of the what the plugin does</source>
        <translation type="obsolete">Zameňte tento text s opisom činnosti daného zásuvného modulu</translation>
    </message>
    <message>
        <source>[menuitemname]</source>
        <translation>[názovpoložkymenu]</translation>
    </message>
    <message>
        <source>&amp;[menuname]</source>
        <translation>&amp;[názovmenu]</translation>
    </message>
    <message>
        <source>Replace this with a short description of what the plugin does</source>
        <translation>Tento text zameňte s krátkym popisom vystihujúcim činnosť vášho zásuvného modulu</translation>
    </message>
</context>
</TS>