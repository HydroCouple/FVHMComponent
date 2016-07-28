#ifndef FVHMCOMPONENT_H
#define FVHMCOMPONENT_H

#include "fvhmcomponent_global.h"
#include "core/headers.h"
#include "core/swmm5.h"
#include "core/abstractmodelcomponent.h"
#include "core/argument1d.h"
#include "core/idbasedargument.h"
#include "temporal/timedata.h"

#include <QDateTime>

class FVHMComponentInfo;
class FVHMInputObjectItem;
class FVHMOutputObjectItem;

class FVHMCOMPONENT_EXPORT FVHMComponent : public AbstractModelComponent
{
      friend class FVHMComponentInfo;

      Q_OBJECT

   public:

      FVHMComponent(const QString &id, FVHMComponent* component = nullptr);

      FVHMComponent(const QString &id, const QString &caption, FVHMComponent* component = nullptr);

      virtual ~FVHMComponent();

      void createArguments();

      void initialize() override;

      void initializeSWMMProject();

      void disposeSWMMProject();

      Project* SWMMProject() const;

      SDKTemporal::Time* startDateTime() const;

      SDKTemporal::Time* endDateTime() const;

      SDKTemporal::Time* currentDateTime() const;

      bool initializeArguments(QString &message);

      void initializeNodeInputExchangeItems();

      void initializeLinkInputExchangeItems();

      void initializeSubCatchmentInputExchangeItems();

      void initializeNodeOutputExchangeItems();

      void initializeLinkOutputExchangeItems();

      void initializeSubCatchmentOutputExchangeItems();

      bool hasError(QString &message);

      HydroCouple::IModelComponent* clone() override;

      QList<QString> validate() override;

      void prepare() override;

      void update(const QList<HydroCouple::IOutput*> &requiredOutputs = QList<HydroCouple::IOutput*>()) override;

      void updateInputExchangeItems();

      void updateOutputExchangeItems();

      void updateOutputExchangeItems(const QList<HydroCouple::IOutput *> &requiredOutputs);

      void finish() override;


   private:
      QList<FVHMInputObjectItem*> m_usedInputs;
      QList<FVHMOutputObjectItem*> m_usedOutputs;
      QStringList m_usedNodes, m_usedLinks, m_usedSubCatchments;
      FVHMComponentInfo  *m_SWMMComponentInfo;
      IdBasedArgumentQString *m_identifiersArgument, *m_inputFilesArgument ;
      Argument1DString *m_usedNodesArgument, *m_usedLinksArgument,*m_usedSubCatchmentsArgument;
      Project* m_SWMMProject;
      QHash<QString,QFileInfo> m_inputFiles;
      SDKTemporal::Time *m_startDateTime, *m_endDateTime, *m_currentDateTime;
      bool m_initialized, m_prepared;
      int m_currentProgress;

};

Q_DECLARE_METATYPE(FVHMComponent*)

#endif // FVHMCOMPONENT_H
